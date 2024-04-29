import argparse
import dash

import pandas as pd
import numpy as np
import plotly.express as px
import dash_bootstrap_components as dbc

app = dash.Dash(external_stylesheets=[dbc.themes.COSMO])

prog_description = "Program used to visualize simulation output from POSE."
parser = argparse.ArgumentParser(prog=prog_description)

parser.add_argument(
    "sim_data", help="Object parameter output generated as output from POSE."
)
parser.add_argument(
    "--object_id",
    default=None,
    required=False,
    help="Object ID to be plotted, not specified if object name is used instead.",
)
parser.add_argument(
    "--object_name",
    default=None,
    required=False,
    help="Object name to be plotted, not specified if object id is used instead.",
)
args = parser.parse_args()

filter_field = None
filter_value = None

if args.object_id is not None:
    filter_field = "id"
    filter_value = int(args.object_id)
elif args.object_name is not None:
    filter_field = "name"
    filter_value = args.object_name
else:
    parser.print_help()
    exit()

sim_data_iter = pd.read_csv(args.sim_data, iterator=True, chunksize=1000)
sim_data = pd.concat(
    [chunk[chunk[filter_field] == filter_value] for chunk in sim_data_iter]
)

sim_time_min = sim_data["sim_time"].min()
sim_time_max = sim_data["sim_time"].max()

app.layout = dash.html.Center(
    children=[
        dash.html.Div(
            children=[
                dash.html.Div(
                    children=[
                        dash.html.Label("Options"),
                        dbc.Checklist(
                            options=[
                                {"label": "Orbit Body", "value": "orbit_body"},
                            ],
                            value=["orbit_body"],
                            id="switch-control",
                            switch=True,
                        ),
                    ],
                    style={"flex": 1, "margin-top": 40},
                ),
                dash.html.Div(
                    children=dash.dcc.Loading(
                        dash.dcc.Graph(
                            id="plot", style={"height": 1000, "width": 1000}
                        ),
                        type="dot",
                    ),
                    style={"flex": 1},
                ),
            ],
            style={"display": "flex", "flexDirection": "row", "width": 1200},
        ),
        dash.html.Div(
            dash.dcc.RangeSlider(
                sim_time_min,
                sim_time_max,
                marks=None,
                value=[sim_time_min, sim_time_max],
                tooltip={
                    "placement": "bottom",
                    "always_visible": True,
                },
                id="sim_time-range",
            ),
            style={"width": 1200},
        ),
    ]
)


def make_sphere(x, y, z, radius, resolution=100):
    """Return the coordinates for plotting a sphere centered at (x,y,z)"""
    u, v = np.mgrid[0 : 2 * np.pi : resolution * 2j, 0 : np.pi : resolution * 1j]
    X = radius * np.cos(u) * np.sin(v) + x
    Y = radius * np.sin(u) * np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)


@dash.callback(
    dash.Output("plot", "figure"),
    dash.Input("switch-control", "value"),
    dash.Input("sim_time-range", "value"),
)
def draw_plot(switch_control, sim_time_range):
    time_mask = (sim_data["sim_time"] >= sim_time_range[0]) & (
        sim_data["sim_time"] <= sim_time_range[1]
    )
    sim_data_slice = sim_data[time_mask]
    orbital_line_fig = px.line_3d(
        sim_data_slice,
        x="x_coord",
        y="y_coord",
        z="z_coord",
        hover_name="name",
        hover_data=["x_coord", "y_coord", "z_coord", "sim_time"],
        title=f"{filter_field}: {filter_value}",
    )

    if "orbit_body" in switch_control:
        orbial_body = make_sphere(0, 0, 0, 6378137.0)
        orbital_line_fig.add_surface(
            x=orbial_body[0],
            y=orbial_body[1],
            z=orbial_body[2],
            hoverinfo="skip",
            colorscale=[[0, "teal"], [1, "teal"]],
            opacity=0.2,
            showscale=False,
        )

    return orbital_line_fig


app.run_server(debug=True)
