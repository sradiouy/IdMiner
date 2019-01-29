import pandas as pd
import itertools
import base64
import datetime
import io
import math
import dash
import dash_table
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
from pathlib import Path
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import networkx as nx
import heapq

from app import app
from src.components.header import headerComponent_termboard

def select_start_query(dfterms):
    try:
        query = dfterms[dfterms["ZIPF_Score"] < 1].sort_values(
            "Genes", ascending=False).iloc[0]["Terms"]
    except:
        query = dfterms.sort_values("Genes", ascending=False).iloc[0]["Terms"]
    return [query]

def get_edges(query,union_intersection_query_result):
    query_name = ",".join(query)
    edge_list = [] #Lista vacia de futuros vertices (termino,gen,numero de articulos donde se menciona el termino y el gen)
    for gene_dict in union_intersection_query_result:
        gene = next(iter(gene_dict))
        articles = gene_dict.get(gene)
        number_articles = float(len(articles.split(","))) # weight del edge
        if number_articles > 0:
            edge_list.append((query_name, gene, number_articles))
    return edge_list

def checknan(value):
    try:
        return math.isnan(float(value))
    except:
        return False

def get_query(term_values,dfgenesbyarticles,union_intersection):
    def get_intersection(genes,terms):
        articles_in_genes = [[str(term[gene]).split(",") for term in terms]for gene in genes]
        list_set = []
        for articles_by_gene in articles_in_genes:
            accumulated_list = set(articles_by_gene[0])
            for s in articles_by_gene[1:]:
                accumulated_list.intersection_update(s)
                list_set.append(accumulated_list)
        return [{gene:",".join(list_set[index])} for index,gene in enumerate(genes) if (len(list_set[index])> 0 and list_set[index] != {"nan"})]
    def get_union(genes,terms):
        dict_gene = {}
        for gene in genes:
            articles = []
            for term in terms:
                if not checknan(term.get(gene)):
                    try:
                        gene_articles = term.get(gene)
                        number_articles = int(len(gene_articles.split(","))) # weight del edge
                    except:
                        gene_articles = str(term.get(gene)).replace(".0","")
                        number_articles = int(len(gene_articles.split(",")))
                    articles.append(gene_articles)
            if len(articles)> 0:
                dict_gene[gene] = ",".join(set(",".join(articles).split(",")))
        return [{key:value} for key,value in dict_gene.items()]
    genes = dfgenesbyarticles.columns[:-1]
    terms  = []
    query = term_values
    if len(query) == 1:
        term = query[0]
        terms = dfgenesbyarticles[dfgenesbyarticles["Terms"] == term].iloc[0][:-1].to_dict()
        return [{key:value} for key,value in terms.items() if not checknan(value)]
    elif len(query) > 1:
        for term in query:
            terms.append(dfgenesbyarticles[dfgenesbyarticles["Terms"] == term].iloc[0][:-1].to_dict())
        if union_intersection == "UNION":
            return get_union(genes,terms)
        elif union_intersection  == "INTERSECTION":
            return get_intersection(genes,terms)
        else:
            raise ValueError("Undefined union_intersection")
    else:
        return False


def get_edges(query,union_intersection_query_result):
        query_name = " - ".join(query)
        edge_list = [] #Lista vacia de futuros vertices (termino,gen,numero de articulos donde se menciona el termino y el gen)
        dict_articles_by_gene = {}
        for gene_dict in union_intersection_query_result:
            gene = next(iter(gene_dict))            
            try:
                articles = gene_dict.get(gene)
                number_articles = int(len(articles.split(","))) # weight del edge
            except:
                articles = str(gene_dict.get(gene)).replace(".0","")
                number_articles = int(len(articles.split(",")))
            if number_articles > 0:
                dict_articles_by_gene[gene] = articles.split(",")
                edge_list.append((query_name, gene, number_articles))
        return edge_list,dict_articles_by_gene


def edge_list_2_article_by_gene(edge_list):
    dict_article_by_gene = {}
    for gene in edge_list:
        dict_article_by_gene[gene[1]] = gene[2] #key: gene-name; value: numb-of-articles
    return dict_article_by_gene



def create_node_trace(network,query,edge_list):
    """Creo un layout del tipo scatter (plotly) para visualizar el network
    
    Returns:
        [go] -- Devuelve un objeto de grafico de plotly
    """
    dict_gene_article = edge_list_2_article_by_gene(edge_list)
    pos = nx.fruchterman_reingold_layout(network)
    node_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers',
        name='ntw',
        hoverinfo='text',
        marker=dict(symbol='circle-dot',
            showscale=True,
            # colorscale options
            # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlOrRd',
            reversescale=True,
            color=[],
            size=[],
            colorbar=dict(
                thickness=30,
                title='Frequency in Articles',
                xanchor='left',
                titleside='right'),
            )
    )
    for node in network.nodes():
        x, y = pos[node]
        node_trace['x'] += tuple([x])
        node_trace['y'] += tuple([y])
    nodecolor = dict([(gene, num_articles["weight"]) for term, gene, num_articles in network.edges(data=True)])
    query_name = " - ".join(query)
    nodecolor[query_name] = 0
    maxinter = max(nodecolor.values())
    for node in nodecolor:
        nodecolor[node] = int((nodecolor[node]/maxinter)*100)
    for node in network.nodes():
        if node != query_name:
            node_info = node + ": " + str(int(dict_gene_article[node])) # Informacion que aparece en el nodo. Numero de articulos por gen.
            node_trace['text'] += tuple([node_info]) # ponemos la informacion en el key text.
            node_trace['marker']['color'] += tuple([int(nodecolor[node])]) # Informacion del escalado. Para colorear por numero de interaccion. El max esta seteado en 100, y es el de mayor interacciones.
            node_trace['marker']['size'] += tuple([70]) #El tamaño del nodo fijo en 70.
        else: # Cuando el nodo es el termino.
            node_trace['marker']['color'] += tuple([0]) # seteo los valores de color 
            node_trace['marker']['size'] += tuple([0]) # y de tamaño a 0
            node_trace['text'] += tuple([""]) # tampoco muesto informacion.
    return node_trace



def create_name_trace(network,query,node_trace,dict_articles_by_gene):
    query_name = " - ".join(query)
    names_nodes = ["<a href=" + "'https://www.ncbi.nlm.nih.gov/pubmed/{0}'".format(",".join(dict_articles_by_gene[gene]).replace(".0","")) + 'style="color: #000000">' + gene + "</a>" if gene != query_name else "" for gene in list(network.nodes())]
    names_trace = go.Scatter(
    x=node_trace["x"],
    y=node_trace["y"],
    text=names_nodes,
    hoverinfo='none',
    textposition='middle center',
    textfont=dict(
        family='arial',
        size=20,
        color='#000000'
    ),
    mode='text')
    return names_trace



def network_layout(query,dict_articles_by_gene):
    query_name = " - ".join(query)
    pubmed = ",".join(set(itertools.chain.from_iterable([value for key,value in dict_articles_by_gene.items()]))).replace(".0","") # Todos los articulos relacionados al )query, en un solo string. Hago el set para que elimine duplicados. Con itertools chains from iterable transofomro una lista de listas a una unica lista
    link = "<a href=" + "'https://www.ncbi.nlm.nih.gov/pubmed/{0}'".format(pubmed) + '>'+ query_name.upper()+'</a>' #Links a todos los articulos
    title = link + ": # Genes: " + str(len(dict_articles_by_gene)) + " ; # Articles: " + str(len(set(pubmed.split(",")))) #Titulo a mostrar
    axis = dict(showline=False,  # hide axis line, grid, ticklabels and  title
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        title=''
        )
    layout = go.Layout(title=title,
        titlefont=dict(
            family='Gadugi',
            size=25,
            color='black'
        ),
        font=dict(size=15),
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
        ),
        annotations=[
            dict(
                showarrow=False,
                text=query_name,
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


def create_network(query,dataframe_genes_by_articles,union_intersection):
    list_articles_by_gene = get_query(query,dataframe_genes_by_articles,union_intersection) # list of dicts, with key gene, value pubmed ids (str concatenated by commas)
    if len(list_articles_by_gene) == 0: #Cuando no hay resultados en una union o interseccion.
        return False
    edge_list,dict_articles_by_gene = get_edges(query,list_articles_by_gene) #Edge-list: list of tuples (term,gene,numb_articles in wich the term appears associated with the gene)
    numb_article_by_gene = edge_list_2_article_by_gene(edge_list) # dict of gene as key and numb_articles in wich the term appears associated with the gene as value.
    network = nx.Graph() # genero el objeto de grafo vacio.
    node_list = [gene[1] for gene in edge_list] # nombre de los genes como nodos
    network.add_weighted_edges_from(edge_list) # Creo lista de arcos. edges.
    node_trace = create_node_trace(network,query,edge_list)
    name_trace = create_name_trace(network,query,node_trace,dict_articles_by_gene)
    layout = network_layout(query,dict_articles_by_gene)
    annot = "<a href='http://www.genomica.weebly.com'>IdMiner: Departamento de Genomica - IIBCE</a>"
    data = [node_trace, name_trace]
    fig = go.Figure(data=data, layout=layout)
    fig['layout']['annotations'][0]['text'] = annot
    return fig


def get_terms_df(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        global dfterms #La seteo global porque como es una varaible que voy a necesitar no hay problemas.
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            dfterms = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')),low_memory=False)
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            dfterms = pd.read_excel(io.BytesIO(decoded),low_memory=False)
        return dfterms
    except Exception as e:
        print(e,filename,"no correct format")


def get_genes_df(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        global dfgenesbyarticles #La seteo global porque como es una varaible que voy a necesitar no hay problemas.
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            dfgenesbyarticles = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            dfgenesbyarticles = pd.read_excel(io.BytesIO(decoded))
        return dfgenesbyarticles
    except Exception as e:
        print(e,filename,"no correct format")
       
    

def parse_contents(dfterms,dfgenesbyarticles):
    try:
        col = dfterms.columns[:-1]
        return html.Div(children=[
            html.Hr(),
            dash_table.DataTable(
            id='table-sorting-filtering',
            columns=[{"name": i, "id": i, 'deletable': True}
                    for i in dfterms[col].columns],
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
                'fontWeight': 'bold',
                'font-family': 'Dosis',
                'textAlign': 'center',
                'font-size': '17'
            },
            style_cell={
                'backgroundColor': '#FAFAFA',
                'minWidth': '0px', 'maxWidth': '150px',
                'whiteSpace': 'no-wrap',
                'overflow': 'hidden',
                'textAlign': 'center',
                'padding': '5px',
                'font-size': '15'
            }),
            html.Hr(),
            dcc.Markdown('''#### Select Query Term:'''),
            html.Div(
                id='term-dropdown-container',
                children=[
                dcc.Dropdown(
                    id='query-term-dropdown',
                    options=[
                        {'label': i.title(), 'value': i} for i in sorted(dfterms.Terms.unique())
                    ],
                    multi=True,
                    value=select_start_query(dfterms) #Start query value
                ),
                dcc.Dropdown(
                    id='union-intersection-dropdown',
                    options=[
                        {'label': 'Union', 'value': 'UNION'},
                        {'label': 'Intersection', 'value': 'INTERSECTION'}
                    ],
                    value='INTERSECTION',
                    searchable=False,
                    clearable=False
                )
                ]
            ),
            html.Hr(),
            dcc.Graph(id='net_graph')
        ])
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])


####main 



uploadOrLoadSample = html.Div(
    className='flex-container',
    children=[
        dcc.Upload(
            id="upload-data",
            className='dashed-file-upload',
            multiple=True,
            children=html.Div(
                ['Drag and Drop or ',html.A('Select Files')]
            )),
        html.Div(id='output-data-upload')
    ]
)

layout = html.Div(
    children=[
        headerComponent_termboard,
        html.Div(
            id="configuration-form-container",
        children=[
            html.H4('EXPLORING TERMS',className='configuration-subsection'),
            uploadOrLoadSample,
            html.Div(id='output-data-upload')
        ]
    )
    ]
)


@app.callback(Output('output-data-upload', 'children'),
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified')])
def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        if len(list_of_contents) ==2:
            for contents, names, dates in zip(list_of_contents, list_of_names, list_of_dates):
                if str(names).split("_")[-1] == "IDMiner-Terms.csv":
                    get_terms_df(contents, names, dates)
                elif names.split("_")[-1] == "IDMiner-Genes.csv":
                    get_genes_df(contents, names, dates)
                else:
                    return html.Div(['There was an error processing this file. You need to put both of this files, (yourname_IDMiner-Terms.csv,yourname_IDMiner-Genes.csv)'])
        else:
             return html.Div(['There was an error processing this file. You need to put BOTH (2) of this files, (yourname_IDMiner-Terms.csv,yourname_IDMiner-Genes.csv)'])
        children = [parse_contents(dfterms, dfgenesbyarticles)]
        return children



@app.callback(
    Output('table-sorting-filtering', 'data'),
    [Input('table-sorting-filtering', 'pagination_settings'),
     Input('table-sorting-filtering', 'sorting_settings'),
     Input('table-sorting-filtering', 'filtering_settings')])
def update_graph(pagination_settings, sorting_settings, filtering_settings):
    filtering_expressions = filtering_settings.split(' && ')
    col = dfterms.columns[:-1]
    dff = dfterms[col]
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


@app.callback(Output('net_graph', 'figure'), [Input('query-term-dropdown', 'value'), 
Input('union-intersection-dropdown', 'value')])
def load_graph(selected_dropdown_value, union_intersection):
    if len(selected_dropdown_value) > 0: #Cuando no tengo seleccion no hago nada. Para evitar que se generen vacios.
        graph = create_network(selected_dropdown_value,dfgenesbyarticles,union_intersection)
        if graph:
            return graph
        else:
            return {} 

