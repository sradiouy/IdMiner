
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from app import app
from src.components import dashboard
from src.components import configuration
from src.components import geneboard



#####


# loads different apps on different urls

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])



@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/':
        return configuration.layout
    elif pathname == '/dashboard':
        return dashboard.layout
    elif pathname == '/geneboard':
        return geneboard.layout
    else:
        return 'The page does not exist: 404 Error'


if __name__ == '__main__':
    app.run_server(debug=False)
