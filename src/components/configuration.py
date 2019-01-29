import os
import time
import logging
import io
import base64

import dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State

from src.components.header import headerComponent_configuration
from src.FetchArticles import *
from src.TermsByAbstracts import *
from src.CleanTerms import *
from src.IdMinerReports import *

from app import app,server


UPLOAD_DIRECTORY = "./project/Results"

if not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)


# Normally, Dash creates its own Flask server internally. By creating our own,
# we can create a route for downloading files directly:


@server.route("/download/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)

def parse_contents(file_name,file_content):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if '.txt' in file_name:
            # Assume that the user uploaded a CSV file
            genefile = open(
                io.StringIO(decoded.decode('utf-8'))).read()
        elif '.fasta' in file_name or '.fa' in file_name or '.faa' in file_name:
            # Assume that the user uploaded an excel file
            genefile = pd.read_excel(io.BytesIO(decoded))
        return genefile
    except Exception as e:
        print(e,file_name,"Uploaded file has not correct format. Extension must be .txt or .fasta!!! \n")


generateNewRun = html.Div(
    children=[
        html.H4('SET UP YOUR RUN', className='configuration-subsection'),
        html.Div(
            className='flex-container',
            children=[
                dcc.Upload(
                    id="upload-file-to-process",
                    className='dashed-file-upload file-upload',
                    children=html.Div(
                        ['Drag and Drop or ',html.A('Select Files')]
                    )
                )
            ],
            title="FASTA or TXT extension!"
        ),
        html.Div(
            className='flex-container term-edition-container',
            children=[
                html.Div(
                    children=[
                        html.Div(
                            children=[
                                html.Label('Keep'),
                                dcc.Textarea(
                                    id='list_keep_terms',
                                    title='Enter (one by line) terms which are frequent in english but you want to keep. See please see documentation to know wich term are in this categroy',
                                    className='area-text',
                                    value='Dominance\nSocial\nLake'
                                )
                            ]
                        ),

                        html.Div(
                            children=[
                                html.Label('Remove'),
                                dcc.Textarea(
                                    id='list_remove_terms',
                                    title='Enter (one by line) terms which are frequent in english but you want to remove. See please see documentation to know wich term are in this categroy',
                                    className='area-text',
                                    value='Human\nPlant\nCancer'
                                ),
                            ]
                        ),
                    ]
                )
            ]
        ),
        html.Div(
            className='box-container',
            children=[
                html.Div(
                    children=[
                        html.Div(
                            children=[
                html.Label('Frequency',title="High frequency means that the word analyzed include word which are frequent in english, like human or cancer.",className="idminer-label"),
                dcc.Dropdown(
                    id='my-term-dropdown',
                    options=[
                        {'label': 'HIGH', 'value': 'high'},
                        {'label': 'MEDIUM', 'value': 'medium'},
                        {'label': 'LOW', 'value': 'low'}
                    ],
                    value="medium",
                    placeholder='Select level...',
                    className="freq-dropdown"
                    )
                        ]
                        ),
                        html.Div(
                            children=[
                                html.Label('Max Terms',className="idminer-label"),
                dcc.Dropdown(
                    id='max-term-dropdown',
                    options=[
                        {'label': '10000', 'value': 10000},
                        {'label': '9000', 'value': 9000},
                        {'label': '8000', 'value': 8000},
                        {'label': '7000', 'value': 7000},
                        {'label': '6000', 'value': 6000},
                        {'label': '5000', 'value': 5000},
                        {'label': '4000', 'value': 4000},
                        {'label': '3000', 'value': 3000},
                        {'label': '2000', 'value': 2000},
                        {'label': '1000', 'value': 1000}
                    ],
                    value=5000,
                    placeholder='Select maximum number of terms...',
                    className="freq-dropdown"
                    ),
                            ]
                        ),
                    ]
                )
            ]
        ),
        html.Div(
            className="box-container",
            children=[
            html.Div(
            children=[
                html.Label('Coverage %',className="idminer-label"),
                dcc.Slider(
                    id='coverage-slider',
                    className='login-slider',
                    min=0,
                    max=100,
                    marks={i: f'{i}%' for i in range(0, 101, 10)},
                    value=50,
                    updatemode='mouseup'
                )
            ]
            ),
            html.Div(
            children=[
                html.Label('Identity %',className="idminer-label"),
                dcc.Slider(
                    id='identity-slider',
                    className='login-slider',
                    min=0,
                    max=100,
                    marks={i: f'{i}%' for i in range(0, 101, 10)},
                    value=50,
                    updatemode='mouseup'
                )
            ]
            )
            ]
        ),
        html.Div(
            className="container4",
            children = [
            html.Button(
                "RUN!",
                id="run-btn",
                title="""
                Please see the wheel spinning at the bottom of the of the page.
                When the run is complete, you would see an status message in the left corner of your screen
                """
            )
            ]
        )
    ]
)

layout = html.Div(
    children=[
        headerComponent_configuration,
        html.Div(
            id='configuration-form-container',
            children=[
            html.Div(
                id='login-container-centered',
                children=[
                    generateNewRun,
                    html.P(id="report-status")
                ]
            ),
            ]
        )
    ]
)







@app.callback(
    Output("report-status", "children"),
    [Input('run-btn', 'n_clicks')],
    [State("upload-file-to-process", "filename"),State("upload-file-to-process", "contents"),State("identity-slider", "value"),State("coverage-slider", "value"),State("max-term-dropdown","value"),State("my-term-dropdown","value"),State("list_keep_terms","value"),State("list_remove_terms","value")],
)
def update_output(clicks,uploaded_filenames, uploaded_file_contents,identity,coverage,maxterms,freqterm,keep,remove):
    """Save uploaded files and regenerate the file list."""
    if "," in keep or "," in remove:
        return [html.P("You must separate keep and remove terms with newline not comma!!!")]
    keep = keep.split("\n") #TERMS TO KEEP TO FILE 
    if len(keep) > 0:
        keep = [term.lower() for term in keep] # TERMS TO KEEP TO FILE 
    remove = remove.split("\n")
    if len(remove) > 0:
        remove = [term.lower() for term in remove]
    print(clicks,uploaded_filenames,identity,coverage,maxterms,freqterm)
    if clicks is not None: #Si cliqueo..
        if uploaded_filenames is None or uploaded_file_contents is None:
            return [html.P("Try to upload the file again!!!")] #No se detecto un archivo
        else:
            if ".fasta" in uploaded_filenames: #if the file has a fasta extension:
                formatfile = "fasta"
                content_type, content_string = uploaded_file_contents.split(',') # get content of the file
                decoded = base64.b64decode(content_string) # decode it  to base64
                input_file = io.StringIO(decoded.decode('utf-8')) # content of the file to strings
            elif ".txt" in uploaded_filenames:
                formatfile = "text"
                content_type, content_string = uploaded_file_contents.split(',') # get content of the file
                decoded = base64.b64decode(content_string) # decode it  to base64
                input_file = io.StringIO(decoded.decode('utf-8'))
            else:
                return [html.P("Format is not .fasta or .txt!!!")]
            if maxterms is None:
                maxterms = 5000
            if freqterm == "high":
                zipf = 5
            elif freqterm == "medium":
                zipf = 3.4
            else:
                zipf = 2 
            if formatfile == "fasta":
                Run_name = "./project/Results/" + uploaded_filenames.split(".fasta")[0]
            else:
                Run_name = "./project/Results/" + uploaded_filenames.split(".txt")[0]
            start_time = time.time()
            logging.basicConfig(format='%(asctime)s - IdMiner - %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',filename= Run_name + '.log',level=logging.DEBUG)
            articles = ids_by_gene(input_file,formatfile,coverage,identity)
            if not articles: #si no obtuvimos ningun articulo
                return [html.P("There were not articles find for your genes. Try to lower de coverage and identity parameters to allow more remote hits")]
            elapsed_time = time.time() - start_time
            logging.info("Fetch article -> Duration %i seconds" %(elapsed_time))
            start_time = time.time()
            gene_pubmed = gene_articles_dict(articles)
            abstractdict = generate_abstracts_dict(gene_pubmed)
            elapsed_time = time.time() - start_time
            logging.info("Abstract by Term -> Duration %i seconds" %(elapsed_time))
            start_time = time.time()
            worddict = cleanwords(abstractdict,keep,remove,zipf)
            dictgram = dict_gram(worddict,maxterms,minfreq=2)
            elapsed_time = time.time() - start_time
            logging.info("Clean Terms -> Duration %i seconds" %(elapsed_time))
            start_time = time.time()
            termdict,geneterms, geneNterms = get_gene_term_dicts(gene_pubmed,worddict,dictgram)
            create_terms_info_dataframe(termdict,Run_name)
            create_gene_artciles_dataframe(geneterms,Run_name)
            elapsed_time = time.time() - start_time
            logging.info("Reports -> Duration %i seconds" %(elapsed_time))
            return [html.P("Done!!!")]
    else:
        return [html.P("")]
    


