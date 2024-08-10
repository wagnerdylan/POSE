import argparse
import json
import numpy
import random
import math
import pandas

from typing import List, Dict

prog_description = "Inject space debris into simulation configuration."
parser = argparse.ArgumentParser(prog=prog_description)
parser.add_argument(
    "sim_config",
    help="Simulation configuration which will be modified in-place with space debris configuration.",
)
parser.add_argument(
    "num_debris",
    help="Number of space debris objects to inject into simulation configuration.",
)
subparsers = parser.add_subparsers(dest='subparser', help="Debris generation modes.")
pure_random_group = subparsers.add_parser("random", description="Pure random debris generation in both coordinate and velocity.")
pure_random_group.add_argument(
    "middle_alt",
    help="Middle altitude to use for the normal distribution of space debris objects.",
)

sim_derived_group = subparsers.add_parser("derived", description="Random debris generation following a trajectory from POSE output.")
sim_derived_group.add_argument("pose_object_file", help="Filepath to a pose object file containing the trajectory for a given object.")
sim_derived_group.add_argument("object_id", help="Object ID of the desired trajectory.")

G = 6.674e-11

EARTH_EQ_RADIUS = 6378137.0
EARTH_MASS = 5.9722e24
DIST_ALT_SD = 1000


def rand_dim(dim: int, scale: float):
    vec = [random.random() for _ in range(dim)]
    mag = sum(x**2 for x in vec) ** 0.5
    unit = [x / mag for x in vec]
    dir = [x * random.choice([1, -1]) for x in unit]
    scaled = [scale * x for x in dir]

    return scaled


def calc_perpendicular_vectors(c_x, c_y, c_z):
    r = math.sqrt(c_x**2 + c_y**2 + c_z**2)
    theta = math.acos(c_z / r)
    phi = math.atan2(c_y, c_x)

    k_x = math.sin(theta) * math.cos(phi)
    k_y = math.sin(theta) * math.sin(phi)
    k_z = math.cos(theta)

    theta += math.pi / 2

    v_x = math.sin(theta) * math.cos(phi)
    v_y = math.sin(theta) * math.sin(phi)
    v_z = math.cos(theta)

    return numpy.asarray([k_x, k_y, k_z]), numpy.asarray([v_x, v_y, v_z])


def rotate_vector_about_k(v, k, theta):
    cos_theta = numpy.cos(theta)
    sin_theta = numpy.sin(theta)

    term1 = v * cos_theta
    term2 = numpy.cross(k, v) * sin_theta
    term3 = k * numpy.dot(k, v) * (1 - cos_theta)

    return term1 + term2 + term3


def tangent_velocity(c_x, c_y, c_z, alt):
    velocity = math.sqrt(G * EARTH_MASS / alt)

    # k pointing up from 0,0,0 toward c_x, c_y, c_z.
    # v is rotated 90deg from k.
    k, v = calc_perpendicular_vectors(c_x, c_y, c_z)

    theta = math.pi / 2
    v_rot = rotate_vector_about_k(v, k, theta)
    v_rot_vel = v_rot * velocity

    return v_rot_vel[0], v_rot_vel[1], v_rot_vel[2]


def sphere_parameters(num: int, radius_cm: float):
    drag_area = math.pi * (radius_cm / 100) ** 2
    drag_coeff = 0.47
    volume = (4 / 3) * math.pi * radius_cm**3
    mass = volume * 2.7  # g/cm^3 density of aluminum
    mass /= 1000
    name = f"d-sph-{num}"

    return name, drag_area, drag_coeff, mass


def create_debris(num: int, c_x: float, c_y: float, c_z: float, v_x: float, v_y: float, v_z: float) -> Dict:
    radius_cm = random.uniform(0.5, 2)
    name, drag_area, drag_coeff, mass = sphere_parameters(num, radius_cm)
    obj = {
        "name": name,
        "drag_area": drag_area,
        "drag_coeff": drag_coeff,
        "mass": mass,
        "radius": radius_cm / 100,
        "state": {
            "soi": {"type": "Earth"},
            "coords": {
                "x": c_x,
                "y": c_y,
                "z": c_z,
            },
            "velocity": {"x": v_x, "y": v_y, "z": v_z},
        },
    }

    return obj


def create_random_debris(num: int, alt: float):
    c_x, c_y, c_z = rand_dim(3, alt)
    v_x, v_y, v_z = tangent_velocity(c_x, c_y, c_z, alt)

    return create_debris(num, c_x, c_y, c_z, v_x, v_y, v_z)


def generate_random_debris(num_debris: int, middle_alt: float) -> List[Dict]:
    assert(middle_alt > 0.0)

    debris = list()
    sc_alt_dist = numpy.random.normal(
        EARTH_EQ_RADIUS + middle_alt, DIST_ALT_SD, num_debris
    )
    for num in range(num_debris):
        debris.append(create_random_debris(num, sc_alt_dist[num]))

    return debris


def rotate_velocity_about_k(c_x, c_y, c_z, v_x, v_y, v_z, theta):
    # Rotated element about k, _v is unused as this term will be calculated using the provided velocity vector.
    k, _v = calc_perpendicular_vectors(c_x, c_y, c_z)
    v_array = numpy.asarray([v_x, v_y, v_z])
    l2_v = numpy.sqrt(numpy.dot(v_array, v_array))
    v = v_array / l2_v

    v_rot = rotate_vector_about_k(v, k, theta)
    return v_rot[0], v_rot[1], v_rot[2]


def create_debris_from_trajectory(num: int, object_trajectory) -> List[Dict]:
    accumulation_factor = num / len(object_trajectory)
    debris = list()
    debris_created = 0
    accumulation_count = 0.0
    
    for point in object_trajectory:
        accumulation_count += accumulation_factor
        num_create = int(accumulation_count // 1)
        for i in range(num_create):
            c_x, c_y, c_z = point["x_coord"], point["y_coord"], point["z_coord"]
            v_x, v_y, v_z = point["x_velocity"], point["y_velocity"], point["z_velocity"]
            
            theta = math.radians(0)
            # 10% chance of a debris object moving in an opposite direction from the defined trajectory.
            if random.random() < 0.1:
                theta = math.radians(180)
            # If more than one object is generated on this point, ensure the trajectory is slightly different.
            elif i > 0:
                theta = math.radians(random.randint(-5,5))
                if theta == 0:
                    theta = 1
            
            vr_x, vr_y, vr_z = rotate_velocity_about_k(c_x, c_y, c_z, v_x, v_y, v_z, theta)
            debris.append(create_debris(debris_created, c_x, c_y, c_z, vr_x, vr_y, vr_z))

            debris_created += 1
            accumulation_count -= 1

    return debris    


def generate_derived_debris(num: int, object_file: str, object_id: int) -> List[Dict]:
    assert(object_id > 0)

    object_df = pandas.read_csv(object_file)
    # Skip over the first part of the object trajectory to prevent collisions on simulation init.
    object_df = object_df[object_df["id"] == object_id]
    object_df = object_df[object_df["sim_time"] > 20.0]
    object_trajectory = object_df.to_dict('records')

    return create_debris_from_trajectory(num, object_trajectory)


def pure_random_generation(args) -> Dict:
    return generate_random_debris(
        int(args.num_debris), float(args.middle_alt)
    )


def sim_derived_generation(args) -> Dict:
    return generate_derived_debris(int(args.num_debris), args.pose_object_file, int(args.object_id))


random.seed(1)
args = parser.parse_args()

with open(args.sim_config, "r") as fd:
    sim_config = json.load(fd)

debris = None
if args.subparser == "random":
    debris = pure_random_generation(args)
elif args.subparser == "derived":
    debris = sim_derived_generation(args)

assert(debris is not None)
sim_config["debris"] = debris

with open(args.sim_config, "w+") as sim_fp:
    json.dump(sim_config, sim_fp, indent=2)
