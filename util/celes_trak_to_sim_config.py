import argparse
import json
import datetime

import dateutil.relativedelta
import pandas as pd

prog_description = "Transform CelesTrak Space Weather Data to a form which may be used by POSE. See https://celestrak.org/SpaceData/."
parser = argparse.ArgumentParser(
    prog="CelesTrak Space Weather Data to POSE Input Configuration"
)
parser.add_argument("input_data", help="CelesTrak Space Weather Data input as csv.")
parser.add_argument(
    "sim_config",
    help="Simulation configuration which to append space weather data to. Start and Halt data must be filled in.",
)
args = parser.parse_args()

input_data = pd.read_csv(
    args.input_data,
    usecols=[
        "DATE",
        "AP1",
        "AP2",
        "AP3",
        "AP4",
        "AP5",
        "AP6",
        "AP7",
        "AP8",
        "F10.7_OBS",
        "F10.7_OBS_CENTER81",
    ],
    converters={"DATE": lambda x: datetime.datetime.fromisoformat(x + "T00:00:00Z")},
)
with open(args.sim_config, "r") as sim_fp:
    sim_config = json.load(sim_fp)

sim_start = datetime.datetime.fromisoformat(sim_config["date"])
sim_stop = datetime.datetime.fromisoformat(sim_config["halt_date"])


def does_date_overlap(d1_start, d1_end, d2_start, d2_end):
    latest_start = max(d1_start, d2_start)
    earliest_end = min(d1_end, d2_end)
    delta = (earliest_end - latest_start).total_seconds() + 1
    overlap = max(0, delta)

    return bool(overlap)


# List of lists, entries are composed of:
#   [0] - Start date of index
#   [1] - AP value for index
#   [2] - F10.7_OBS value for index
#   [3] - F10.7_OBS_CENTER81 value for index
data_list = []
for index, row in input_data.iterrows():
    if pd.isnull(row["AP1"]):
        interval_end = row["DATE"] + dateutil.relativedelta.relativedelta(months=+1)
        if does_date_overlap(
            row["DATE"],
            interval_end,
            sim_start,
            sim_stop,
        ):
            data_list.append(
                [
                    row["DATE"].isoformat(),
                    interval_end.isoformat(),
                    5.0,
                    row["F10.7_OBS"],
                    row["F10.7_OBS_CENTER81"],
                ]
            )
    else:
        for ap_i in range(8):
            interval_start = row["DATE"] + datetime.timedelta(hours=3 * ap_i)
            interval_end = interval_start + dateutil.relativedelta.relativedelta(hours=+3)
            if does_date_overlap(
                interval_start,
                interval_end,
                sim_start,
                sim_stop,
            ):
                data_list.append(
                    [
                        interval_start.isoformat(),
                        interval_end.isoformat(),
                        row["AP" + str(ap_i + 1)],
                        row["F10.7_OBS"],
                        row["F10.7_OBS_CENTER81"],
                    ]
                )

# Update simulation configuration with space weather data for earth.
sim_config["earth_sw"] = data_list

# Write-out new simulation configuration in-place.
with open(args.sim_config, "w+") as sim_fp:
    json.dump(sim_config, sim_fp, indent=2)
