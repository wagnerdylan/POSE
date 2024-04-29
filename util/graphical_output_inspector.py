import argparse
import dash

import plotly.graph_objects as pgo
import pandas as pd
import numpy as np
import dash_bootstrap_components as dbc

app = dash.Dash(external_stylesheets=[dbc.themes.COSMO])

prog_description = "Program used to visualize simulation output from POSE."
parser = argparse.ArgumentParser(prog=prog_description)

parser.add_argument(
    "sim_data", help="Object parameter output generated as output from POSE."
)
parser.add_argument(
    "--object_id",
    default=[],
    required=False,
    action="append",
    help="Object ID to be plotted, option may be repeated to plot multiple objects.",
)
parser.add_argument(
    "--object_name",
    default=[],
    required=False,
    action="append",
    help="Object name to be plotted, option may be repeated to plot multiple objects.",
)
args = parser.parse_args()
plot_objects = [{"filter": "name", "value": x} for x in args.object_name] + [
    {"filter": "id", "value": int(x)} for x in args.object_id
]

if not len(plot_objects):
    raise ValueError("At least one object must be specified as input.")

sim_data = pd.read_csv(args.sim_data)

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


def format_text_info(row):
    hover_info = [
        "name",
        "sim_time",
        "x_coord",
        "y_coord",
        "z_coord",
        "altitude",
        "x_velocity",
        "y_velocity",
        "z_velocity",
    ]
    text = ""
    for attr in hover_info:
        text += f"{attr}: {row[attr]}<br>"
    return text


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

    orbital_line_fig = pgo.Figure()

    for plot_object in plot_objects:
        object_data = sim_data_slice[
            sim_data_slice[plot_object["filter"]] == plot_object["value"]
        ]
        text_labels = [format_text_info(row) for _, row in object_data.iterrows()]
        orbital_line_fig.add_trace(
            pgo.Scatter3d(
                x=object_data["x_coord"],
                y=object_data["y_coord"],
                z=object_data["z_coord"],
                text=text_labels,
                hoverinfo="text",
                name=object_data["name"].values[0],
                mode="lines",
            )
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
