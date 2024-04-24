import argparse
import json
import numpy
import random
import math

from typing import List, Dict

prog_description = "Inject space debris into simulation configuration."
parser = argparse.ArgumentParser(prog=prog_description)
parser.add_argument("sim_config", help="Simulation configuration which will be modified in-place with space debris configuration.")
parser.add_argument("num_debris", help="Number of space debris objects to inject into simulation configuration.")
parser.add_argument("middle_alt", help="Middle altitude to use for the normal distribution of space debris objects.")

args = parser.parse_args()

G = 6.674e-11

EARTH_EQ_RADIUS = 6378137.0
EARTH_MASS = 5.9722e24
DIST_ALT_SD = 1000

def rand_dim(dim: int, scale: float):
    vec = [random.random() for _ in range(dim)]
    mag = sum(x**2 for x in vec) ** .5
    unit = [x/mag for x in vec]
    dir = [x*random.choice([1, -1]) for x in unit]
    scaled = [scale*x for x in dir]

    return scaled

def tangent_velocity(c_x, c_y, c_z, alt):
    velocity = math.sqrt(G*EARTH_MASS / alt)

    r = math.sqrt(c_x**2 + c_y**2 + c_z**2)
    theta = math.acos(c_z/r)
    phi = math.atan2(c_y, c_x)

    theta += math.pi / 2

    v_x = math.sin(theta) * math.cos(phi) * velocity
    v_y = math.sin(theta) * math.sin(phi) * velocity
    v_z = math.cos(theta) * velocity

    return v_x, v_y, v_z

def sphere_parameters(num: int, radius_cm: float):
    drag_area = math.pi * (radius_cm/100)**2
    drag_coeff = 0.47
    volume = (4/3) * math.pi * radius_cm**3
    mass = volume * 2.7 # g/cm^3 density of aluminum
    mass /= 1000
    name = f"d-sph-{num}"

    return name, drag_area, drag_coeff, mass

def random_debris(num: int, alt: float) -> Dict:
    c_x, c_y, c_z = rand_dim(3, alt)
    v_x, v_y, v_z = tangent_velocity(c_x, c_y, c_z, alt)
    name, drag_area, drag_coeff, mass = sphere_parameters(num, random.uniform(0.5, 2))
    obj = {
        "name": name,
        "soi": {
            "type": "Earth"
        },
        "coords": {
            "x": c_x,
            "y": c_y,
            "z": c_z,
        },
        "velocity": {
            "x": v_x,
            "y": v_y,
            "z": v_z
        },
        "drag_area": drag_area,
        "drag_coeff": drag_coeff,
        "mass": mass
    }

    return obj

def generate_random_debris(num_debris: int, middle_alt: float) -> List[Dict]:
    debris = list()
    sc_alt_dist = numpy.random.normal(EARTH_EQ_RADIUS + middle_alt, DIST_ALT_SD, num_debris)
    for num in range(num_debris):
        debris.append(random_debris(num, sc_alt_dist[num]))
    
    return debris

with open(args.sim_config, "r") as fd:
    sim_config = json.load(fd)

sim_config["debris"] = generate_random_debris(int(args.num_debris), int(args.middle_alt))

with open(args.sim_config, "w+") as sim_fp:
    json.dump(sim_config, sim_fp, indent=2)