import plotly.plotly as py
import plotly.graph_objs as go

import plotly.tools as tls
tls.set_credentials_file(username='annevanrossum', api_key='2z29vr6xo7')

import numpy

import pandas as pd

# Read data from a csv
#z_data = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/api_docs/mt_bruno_elevation.csv')
# Treat first row as header names
z_data = pd.read_csv('angular_line.plotly.csv', header=0, index_col=0)

data = [
    go.Surface(
        x=z_data.axes[0],
        y=z_data.axes[1],
        z=z_data.as_matrix(),
        showscale=False
    )
]
layout = go.Layout(
    title='Angular Line',
    autosize=False,
    width=500,
    height=500,
    margin=dict(
        l=65,
        r=50,
        b=65,
        t=90
    )
)
fig = go.Figure(data=data, layout=layout)
py.plot(fig, filename='angular-line')

