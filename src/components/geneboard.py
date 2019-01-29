#GeneBoard. In common-gene 2 gene. 
import base64
import datetime
import io
import os 
import pandas as pd 
import itertools
from collections import defaultdict
import networkx as nx 

import plotly.graph_objs as go
import dash
import dash_table
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html

from app import app
from src.components.header import headerComponent_geneboard


def get_genes_relations(dfgenes):
    """ Esta funcion busca obtener a partir de la table de articulos por termino por genes, la cantidad de articulos asociados a un gen particular, y la cantidad de articulos comunes entre los genes.
    
    Arguments:
        dfgenes {[dataframe]} -- dfgenes: es un dataframe resultado de la corrida de IDMINER.
    
    Returns:
        [dict] -- [article_by_gene: diccionario donde el key es el gen y el value es la lista de articulos relacionados al gen]
        [dict] -- [gene_common: diccionario donde el key es el gen y el value es la lista de diccionario donde los key son los genes (distintos al primero) y los values son los articulos que tienen los genes en comun]
        [dict] -- [gene_common_count: diccionario donde el key es el gen y el value es la lista de diccionario donde los key son los genes (distintos al primero) y los values son la cantidad de articulos que tienen los genes en comun]
    """

    genes = dfgenes.columns[0:-1]
    article_by_gene = {}
    for gene in genes:
        article_by_gene[gene] = list(set(itertools.chain.from_iterable(dfgenes[gene].dropna().astype(str).str.split(",").tolist())))

    gene_common = defaultdict(dict)
    gene_common_count = defaultdict(dict)

    for gene in genes:
        for gene2 in genes:
            if gene != gene2:
                in_common = list(set(article_by_gene[gene]) & set(article_by_gene[gene2]))
                if len(in_common) > 0:
                    gene_common[gene][gene2] = in_common
                    gene_common_count[gene][gene2] = len(in_common)

    return article_by_gene,gene_common,gene_common_count



def generate_dataframe(article_by_gene,gene_common_count):
    """[summary]
    
    Arguments:
        article_by_gene {[dict]} -- Diccionario donde el key es el gen y el value es la lista de articulos relacionados al gen
        gene_common_count {[dict]} -- Diccionario donde el key es el gen y el value es la lista de diccionario donde los key son los genes (distintos al primero) y los values son la cantidad de articulos que tienen los genes en comun
    
    Returns:
        [dataframe] -- Dataframe que muestra el numero de articulos en comun entre dos genes.
    """

    genes = list(article_by_gene.keys()) #Defino los nombres de los genes para que sean las futuras columnas
    rows = []
    start_query_value = 0
    start_query_max = 0
    start_query = ""
    for first_gene in genes: #Obtengo las relaciones entre los genes
        row = []
        for second_gene in genes:
            if first_gene == second_gene:
                own_articles =len(article_by_gene[first_gene])
                row.append(own_articles)
            else:
                if second_gene in gene_common_count[first_gene].keys():
                    row.append(gene_common_count[first_gene][second_gene])
                else:
                    row.append(0)
        start_query_value = sum(row) - own_articles
        if start_query_value > start_query_max:
            start_query_max = start_query_value 
            start_query = first_gene
        row.insert(0,len(gene_common_count[first_gene]))
        row.insert(0, first_gene)            
        rows.append(row)
    genes.insert(0,"In-common")    
    genes.insert(0,"Gene")
    df = pd.DataFrame(rows, columns=genes)
    return df,start_query


def get_gene_edges(query_gene,gene_common_count):
    """Funcion que me permite crear la estrucutra de tuplas (lista de tuplas), necesario para establecer los edges del grafo. Cada tupla consiste en el query_gene, el subject_gene y la cantidad de articulos comunes a estos dos. 
    
    Arguments:
        query_gene {[str]} -- Nombre del gen del cual quiero establecer relaciones
        gene_common_count {[dict]} -- diccionario donde el key es el gen y el value es la lista de diccionario donde los key son los genes (distintos al primero) y los values son la cantidad de articulos que tienen los genes en comun
    
    Returns:
        [list] -- [Lista de tuplas de 3 valores: Query_gene (nodo 1); Subject_gene (nodo 2) y numero de articulos en comun (peso de la interaccion)]
    """

    gene_dict = gene_common_count[query_gene] #Obtengo el diccionario para el gen de interes
    edge_list = [] #Lista vacia de futuros vertices (termino,gen,numero de articulos donde se menciona el termino y el gen)
    for subject_gene,num_articles in gene_dict.items(): #itero el diccionario para establecer las relaciones con el resto (subjects) de los genes.
        edge_list.append((query_gene, subject_gene, num_articles))
    return edge_list


def create_gene_node_trace(network,query_gene,edge_list,gene_common_count):
    """Creo un layout del tipo scatter (plotly) para visualizar el network. En este caso defino el lugar donde van a estar los nodos en un objeto del tipo scatter.
    
    Returns:
        [go] -- Devuelve un objeto de grafico scatter de plotly. 
    """
    pos = nx.fruchterman_reingold_layout(network) #Defino un tipo de layout, en este caso elijo el reingold.
    node_trace = go.Scatter(
        x=[], #posicion en el eje x del nodo
        y=[],#posicion en el eje y del nodo
        text=[],
        mode='markers', #marcadores 
        name='ntw',
        hoverinfo='text', # informacion que me da al pasar por arriba del punto
        marker=dict(symbol='circle-dot',
            showscale=True, # que me muestre escala 
            colorscale='YlOrRd', # escala de colores: amarillo (low) al rojo (high)
            reversescale=True, 
            color=[], # Lista donde va a estar codificado los colores de cada nodo
            size=[], # Lista donde va a estar codificado el tamaño de cada nodo
            colorbar=dict( # barra de color para mostrar la escala de color. 
                thickness=30,
                title='Articles in common', 
                xanchor='left',
                titleside='right'),
            )
    )
    for node in network.nodes(): #itero los nodos de mi grafo
        x, y = pos[node] # obtengo las posiciones a partir del mapa generado por reingold
        node_trace['x'] += tuple([x]) # posicion en el eje x del nodo
        node_trace['y'] += tuple([y]) # posicion en el eje y del nodo
    nodecolor = dict([(gene, num_articles["weight"]) for term, gene, num_articles in network.edges(data=True)]) #Creo un diccionario para determinar como indicador del color el numero de genes en comun entre los gene subjects y el gene query.
    nodecolor[query_gene] = 0 # determino que el gene-query no tenga color (no quiero que aparezca en el grafo)
    maxinter = max(nodecolor.values()) # determino como valor maximo de escala el gene_query:gene_subject que tenga mas articulos en comun
    for node in nodecolor: #normalizo por ese valor
        nodecolor[node] = int((nodecolor[node]/maxinter)*100)
    for node in network.nodes(): #iterno los nodos del network (grafo)
        if node != query_gene: # en el caso que sea distinto al query-gene
            node_info = node + ": " + str(gene_common_count[query_gene][node]) # Informacion que aparece en el nodo. Numero de articulos por gen.
            node_trace['text'] += tuple([node_info]) # ponemos la informacion en el key text.
            node_trace['marker']['color'] += tuple([int(nodecolor[node])]) # Informacion del escalado. Para colorear por numero de interaccion. El max esta seteado en 100, y es el de mayor interacciones.
            node_trace['marker']['size'] += tuple([70]) #El tamaño del nodo fijo en 70.
        else: # Cuando el nodo es el gene-query.
            node_trace['marker']['color'] += tuple([0]) # seteo los valores de color 
            node_trace['marker']['size'] += tuple([0]) # y de tamaño a 0
            node_trace['text'] += tuple([""]) # tampoco muesto informacion.
    return node_trace




def create_gene_name_trace(network,query_gene,node_trace,gene_common):
    """[summary]
    
    Arguments:
        network {[networx object]} -- Grafo creado a partir de networkx 
        query_gene {[str]} -- Gen de interes
        node_trace {[dict]} -- Diccionario que contiene informacion de los nodos en un objeto del tipo scatter de plotly.
        gene_common {[dict]} -- Diccionario donde el key es el gen y el value es la lista de diccionario donde los key son los genes (distintos al primero) y los values son los articulos que tienen los genes en comun
    
    Returns:
        [dict] -- [names_nodes: diccionario que contiene la informacion del nombre de los nodos y articulos que se asocian a cada nodo(gen)]
    """

    names_nodes = ["<a href=" + "'https://www.ncbi.nlm.nih.gov/pubmed/{0}'".format(",".join(gene_common[query_gene][gene]).replace(".0","")) + 'style="color: #000000">' + gene + "</a>" if gene != query_gene else "" for gene in list(network.nodes())] #Determino el nombre de los nodos. El nombre de los nodos esta compuesto por el nombre del gene subject y contiene un link a los articulos que tiene en comun con el gene-query.
    names_trace = go.Scatter( #creo un objeto scatter
    x=node_trace["x"], #defino la posicion de los nombres en el lugar donde previamente determino los nodos (node_trace)
    y=node_trace["y"], #defino la posicion de los nombres en el lugar donde previamente determino los nodos (node_trace)
    text=names_nodes, #lista que contiene el nombre de los nodos
    hoverinfo='none', # no agregar informacion del tipo hover (no pasa nada si paso por arriba), ya que previamente determine esa informacion en el node_trace
    textposition='middle center', # lo ubico en el centro.
    textfont=dict( #determino el tipo de letra.
        family='arial',
        size=20,
        color='#000000'
    ),
    mode='text')
    return names_trace


def gene_network_layout(query_gene,gene_common,articles_by_gene):
    """Genero el layout del grafo.
    
    Arguments:
        query_gene {[str]} --  Nombre del gen del cual quiero establecer relaciones
        articles_by_gene {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """

    pubmed = ",".join(article_by_gene[query_gene]).replace(".0","")  # Todos los articulos relacionados al query, en un solo string. Hago el set para que elimine duplicados. Con itertools chains from iterable transofomro una lista de listas a una unica lista
    link = "<a href=" + "'https://www.ncbi.nlm.nih.gov/pubmed/{0}'".format(pubmed) + '>'+ query_gene +'</a>' #Links a todos los articulos
    title = link + ": # Genes: " + str(len(gene_common[query_gene]))  + "; # Articles: " + str(len(article_by_gene[query_gene])) #Determino titulo
    axis = dict(showline=False,  # hide axis line, grid, ticklabels and  title
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        title=''
        )
    layout = go.Layout(title=title, # Establezco estilo del titulo
        titlefont=dict(
            family='Gadugi',
            size=25,
            color='black'
        ),
        font=dict(size=15), #Establezco estilo del grafo.
        plot_bgcolor='#EDEEF0',
        showlegend=False,
        autosize=True,
        height=800,
        xaxis=go.layout.XAxis(axis),
        yaxis=go.layout.YAxis(axis),
        margin=go.layout.Margin(
            l=100,
            r=100,
            b=100,
            t=100,
        ), #determino estilo de anotacio
        annotations=[
            dict(
                showarrow=False,
                text=query_gene,
                xref='paper',
                yref='paper',
                x=0,
                y=-0.1,
                xanchor='left',
                yanchor='bottom',
                font=dict(
                    size=20
                )
            )
        ]
        )
    return layout


def create_gene_network(query_gene,article_by_gene,gene_common,gene_common_count):
    """Creacion del grafo. Centrado en un gen y la realcion con el resto de los genes (articulos en comun)
    
    Arguments:
        query_gene {[str]} --  Nombre del gen del cual quiero establecer relaciones
        article_by_gene {[dict]} -- Diccionario donde el key es el gen y el value es la lista de articulos relacionados al gen
        gene_common {[dict]} -- Diccionario donde el key es el gen y el value es la lista de diccionario donde los key son los genes (distintos al primero) y los values son los articulos que tienen los genes en comun
        gene_common_count {[dict]} -- Diccionario donde el key es el gen y el value es la lista de diccionario donde los key son los genes (distintos al primero) y los values son la cantidad de articulos que tienen los genes en comun
    
    Returns:
        [dict] -- fig: es un diccionario que contiene la figura del grafo, que es ploteado por dassh-plotly
    """

    edge_list = get_gene_edges(query_gene,gene_common_count) # Creo la estructura de EDGES.
    if len(edge_list) == 0: #Cuando no hay resultados en una union o interseccion.
        return False
    network = nx.Graph() # genero el objeto de grafo vacio.
    node_list = [gene[1] for gene in edge_list] # nombre de los genes como nodos
    network.add_weighted_edges_from(edge_list) # Creo lista de arcos. edges.
    node_trace = create_gene_node_trace(network,query_gene,edge_list,gene_common_count)
    name_trace = create_gene_name_trace(network,query_gene,node_trace,gene_common)
    layout = gene_network_layout(query_gene,gene_common,article_by_gene)
    annot = "<a href='http://www.genomica.weebly.com'>IdMiner: Departamento de Genomica - IIBCE</a>"
    data = [node_trace, name_trace]
    fig = go.Figure(data=data, layout=layout)
    fig['layout']['annotations'][0]['text'] = annot
    return fig


def get_genes_df(contents, filename, date):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    global df_genes_relation
    global start_gene_query
    global article_by_gene
    global gene_common
    global gene_common_count
    global dfgenes
    try:
         #La seteo global porque como es una varaible que voy a necesitar no hay problemas.
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            dfgenes = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            dfgenes = pd.read_excel(io.BytesIO(decoded))
        article_by_gene,gene_common,gene_common_count = get_genes_relations(dfgenes)
        df_genes_relation, start_gene_query = generate_dataframe(article_by_gene,gene_common_count)
        return None
    except Exception as e:
        print(e,filename,"no correct format")



def parse_contents(df_genes_realtion):
    try:
        col = df_genes_realtion.columns
        return html.Div(children=[
            html.Hr(),
            dash_table.DataTable(
                id='gene_table-sorting-filtering',
                columns=[{"name": i, "id": i, 'deletable': True}
                        for i in df_genes_realtion[col].columns],
                pagination_settings={
                    'current_page': 0,
                    'page_size': 10
                },
                pagination_mode='be',
                sorting='be',
                sorting_type='multi',
                sorting_settings=[],
                filtering='be',
                filtering_settings='',
                style_table={'overflowX': 'scroll'},
                style_header={
                    'backgroundColor': '#91B9E5',
                    'minWidth': '0px', 'maxWidth': '800px',
                    'fontWeight': 'bold',
                    'font-family': 'inherit',
                    'textAlign': 'center',
                    'padding-right': '20px',
                    'padding-left': '20px',
                    'font-size': '15'
                },
                style_cell={
                    'backgroundColor': '#FAFAFA',
                    'minWidth': '0px', 'maxWidth': '800px',
                    'whiteSpace': 'no-wrap',
                    'overflow': 'hidden',
                    'textAlign': 'center',
                    'font-size': '15'
                }
            ),
            html.Hr(),
            dcc.Markdown('''#### Select Gene:'''),
            html.Div(
                id='gene-dropdown-container',
                children=[
                    dcc.Dropdown(
                        id='query-gene-dropdown',
                        options=[
                            {'label': i.title(), 'value': i} for i in sorted(dfgenes.columns[:-1])
                        ],
                        multi=False,
                        value=start_gene_query #Start query value
                    )
                ]
            ),
        html.Hr(),
        dcc.Graph(id='gene_net_graph')
        ]
        )
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])





#MAIN



uploadOrLoadSample = html.Div(
    className='flex-container',
    children=[
        dcc.Upload(
            id="upload-gene-data",
            className='dashed-file-upload',
            children=html.Div(
                ['Drag and Drop or ',html.A('Select Files')]
            )),
        html.Div(id='output-gene-data-upload')
    ]
)


layout = html.Div(
    children=[
        headerComponent_geneboard,
        html.Div(
            id="configuration-form-container",
        children=[
            html.H4('EXPLORING GENES',className='configuration-subsection'),
            uploadOrLoadSample,
            html.Div(id='output-gene-data-upload')
        ]
    )
    ]
)


@app.callback(Output('output-gene-data-upload', 'children'),
              [Input('upload-gene-data', 'contents')],
              [State('upload-gene-data', 'filename'),
               State('upload-gene-data', 'last_modified')])
def update_output(file_content, file_name, file_date):
    if file_content is not None:
        if file_name.split("_")[-1] == "IDMiner-Genes.csv":
            get_genes_df(file_content, file_name, file_date)
            children = [parse_contents(df_genes_relation)]
        else:
             return html.Div(['There was an error processing this file. You need to upload this file, (yourname_IDMiner-Genes.csv)'])
        return children


@app.callback(
    Output('gene_table-sorting-filtering', 'data'),
    [Input('gene_table-sorting-filtering', 'pagination_settings'),
     Input('gene_table-sorting-filtering', 'sorting_settings'),
     Input('gene_table-sorting-filtering', 'filtering_settings')])
def update_graph(pagination_settings, sorting_settings, filtering_settings):
    filtering_expressions = filtering_settings.split(' && ')
    dff = df_genes_relation
    for filter in filtering_expressions:
        if ' eq ' in filter:
            col_name = filter.split(' eq ')[0]
            filter_value = filter.split(' eq ')[1]
            dff = dff.loc[dff[col_name] == filter_value]
        if ' > ' in filter:
            col_name = filter.split(' > ')[0]
            filter_value = float(filter.split(' > ')[1])
            dff = dff.loc[dff[col_name] > filter_value]
        if ' < ' in filter:
            col_name = filter.split(' < ')[0]
            filter_value = float(filter.split(' < ')[1])
            dff = dff.loc[dff[col_name] < filter_value]

    if len(sorting_settings):
        dff = dff.sort_values(
            [col['column_id'] for col in sorting_settings],
            ascending=[
                col['direction'] == 'asc'
                for col in sorting_settings
            ],
            inplace=False
        )

    return dff.iloc[
        pagination_settings['current_page']*pagination_settings['page_size']:
        (pagination_settings['current_page'] + 1) *
        pagination_settings['page_size']
    ].to_dict('rows')


@app.callback(Output('gene_net_graph', 'figure'), [Input('query-gene-dropdown', 'value')])
def load_graph(selected_dropdown_value):
    if len(selected_dropdown_value) > 0: #Cuando no tengo seleccion no hago nada. Para evitar que se generen vacios.
        graph = create_gene_network(selected_dropdown_value,article_by_gene,gene_common,gene_common_count)
        if graph:
            return graph
        else:
            return {} 