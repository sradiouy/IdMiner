import dash_html_components as html
import dash_core_components as dcc



def get_menu():
    menu = html.Div(children=[

        dcc.Link('Term Exploration    ', href='/dashboard', className="tab first"),
        html.Br(),
        dcc.Link('Gene Exploration    ', href='/geneboard', className="tab")

    ], className="row inline-flex-container")
    return menu

def get_menu_term():
    menu = html.Div(children=[

        dcc.Link('Configuration    ', href='/', className="tab first"),
        html.Br(),
        dcc.Link('Gene Exploration    ', href='/geneboard', className="tab")

    ], className="row inline-flex-container")
    return menu


def get_menu_gene():
    menu = html.Div(children=[

        dcc.Link('Configuration    ', href='/', className="tab first"),
        html.Br(),
        dcc.Link('Term Exploration    ', href='/dashboard', className="tab")

    ], className="row inline-flex-container")
    return menu

headerComponent_configuration = html.Div(
    children=[
        html.Div(
            children=[
            html.Img(src='https://raw.githubusercontent.com/sradiouy/IdMiner/master/logo_transparent_background.png',
                style={
                    'height': '250px',
                    "margin-left": "auto",
	                "margin-right": "auto",
	                "display": "block",
                },
            ),
            ], 
        className='header-container'
        ),
        html.Br([]),
        get_menu()
        ]
        )


headerComponent_termboard = html.Div(
    children=[
        html.Div(
            children=[
            html.Img(src='https://raw.githubusercontent.com/sradiouy/IdMiner/master/logo_transparent_background.png',
                style={
                    'height': '200px',
                    "margin-left": "auto",
	                "margin-right": "auto",
	                "display": "block",
                },
            ),
            ], 
        className='header-container'
        ),
        html.Br([]),
        get_menu_term()
        ]
        )

headerComponent_geneboard = html.Div(
    children=[
        html.Div(
            children=[
            html.Img(src='https://raw.githubusercontent.com/sradiouy/IdMiner/master/logo_transparent_background.png',
                style={
                    'height': '200px',
                    "margin-left": "auto",
	                "margin-right": "auto",
	                "display": "block",
                },
            ),
            ], 
        className='header-container'
        ),
        html.Br([]),
        get_menu_gene()
        ]
        )



