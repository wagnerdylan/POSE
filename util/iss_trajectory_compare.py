import argparse
import pandas as pd
import numpy as np
import datetime

parser = argparse.ArgumentParser(prog="ISS Trajectory Compare")
parser.add_argument(
    "iss_data",
    help="Filepath to the iss trajectory file to be used for comparison.",
)
parser.add_argument(
    "sim_data",
    help="Filepath to the simulation output file for objects.",
)
parser.add_argument("output", help="Destination to write resulting table diff to.")
args = parser.parse_args()

sim_data = pd.read_csv(
    args.sim_data,
    usecols=[
        "sim_time",
        "x_coord",
        "y_coord",
        "z_coord",
        "x_velocity",
        "y_velocity",
        "z_velocity",
    ],
    converters={
        "sim_time": lambda x: (
            datetime.datetime.fromisoformat("2023-10-23T12:00:00.000")
            + datetime.timedelta(seconds=float(x))
        )
    },
)
conv_to_meters = lambda x: np.float64(x) * 1000
iss_data = pd.read_csv(
    args.iss_data,
    sep=" ",
    skiprows=38,
    header=None,
    names=[
        "sim_time",
        "x_coord",
        "y_coord",
        "z_coord",
        "x_velocity",
        "y_velocity",
        "z_velocity",
    ],
    converters={
        "sim_time": lambda x: datetime.datetime.fromisoformat(x),
        "x_coord": conv_to_meters,
        "y_coord": conv_to_meters,
        "z_coord": conv_to_meters,
        "x_velocity": conv_to_meters,
        "y_velocity": conv_to_meters,
        "z_velocity": conv_to_meters,
    },
)

sim_data = sim_data[sim_data["sim_time"].isin(iss_data["sim_time"])].reset_index(
    drop=True
)
diff = sim_data.subtract(iss_data)
diff["sim_time"] = sim_data["sim_time"]
diff = diff.dropna()
diff.to_csv(args.output)
