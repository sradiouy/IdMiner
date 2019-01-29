import dash
from flask import Flask, send_from_directory


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']



# Normally, Dash creates its own Flask server internally. By creating our own,
# we can create a route for downloading files directly:
server = Flask(__name__)

app = dash.Dash(__name__, external_stylesheets=external_stylesheets,server=server)
app.config.supress_callback_exceptions = True
