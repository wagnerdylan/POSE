import json
import random
import sys
import getopt
import uuid
import datetime
import math

# Better Spherical Random Point Function

radus_of_earth = 6378100  # in meters


def get_point():
    u = random.random()
    v = random.random()
    theta = u * 2.0 * math.pi
    phi = math.acos(2.0 * v - 1.0)
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)
    x = radus_of_earth * sin_phi * cos_theta
    y = radus_of_earth * sin_phi * sin_theta
    z = radus_of_earth * cos_phi
    return x, y, z


# Define Low Orbit Statistics(Circular Orbit) https://en.wikipedia.org/wiki/Orbital_speed

min_low_orbit = 6578100.0  # m 200,000 min w/o radius
max_low_orbit = 8378100.0  # m 2,000,000 max w/o radius

minimal_velocity = 6900.0  # m/s
max_velocity = 7800.0  # m/s

# Linear Scale, might want to switch to an exponential scale


def scaling_velocity(x: float) -> float:
    return (((max_velocity-minimal_velocity)*(x-min_low_orbit))/(max_low_orbit - min_low_orbit)) + minimal_velocity


data = {"date": "", "debris": [], "spacecraft": []}
ref = "data/"

"""make distance respect earth's radius, make velocity believable"""


def rand_orbit() -> float:
    return random.uniform(min_low_orbit, max_low_orbit)


def create_spacecraft() -> dict:
    spacecraft = {}
    # spacecraft["id"] = uuid.uuid4().int & (1 >> 32)-1
    spacecraft["id"] = random.getrandbits(32)
    rand_orbital = rand_orbit()
    rand_x, rand_y, rand_z = get_point()
    spacecraft["x_dis"] = rand_x + rand_orbital
    spacecraft["y_dis"] = rand_y + rand_orbital
    spacecraft["z_dis"] = rand_z + rand_orbital
    spacecraft["x_vel"] = scaling_velocity(rand_x + rand_orbital)
    spacecraft["y_vel"] = scaling_velocity(rand_y + rand_orbital)
    spacecraft["z_vel"] = scaling_velocity(rand_z + rand_orbital)
    return spacecraft


def create_debris() -> dict:
    debris = {}
    # debris["id"] = uuid.uuid4().int & (1 >> 32)-1
    debris["id"] = random.getrandbits(32)
    rand_orbital = rand_orbit()
    rand_x, rand_y, rand_z = get_point()
    debris["x_dis"] = rand_x + rand_orbital
    debris["y_dis"] = rand_y + rand_orbital
    debris["z_dis"] = rand_z + rand_orbital
    debris["x_vel"] = scaling_velocity(rand_x + rand_orbital)
    debris["y_vel"] = scaling_velocity(rand_y + rand_orbital)
    debris["z_vel"] = scaling_velocity(rand_z + rand_orbital)
    return debris


def main(argv):

    output_file = ''
    debris_num = 0
    spacecraft_num = 0

    try:
        opts, args = getopt.getopt(argv, "hd:s:o:")
    except getopt.GetoptError:
        print('generate.py -h <help> -d <number of debris> -s <number of spacecraft> -o <output file name>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('generate.py -h <help> -d <number of debris> -s <number of spacecraft> -o <output file name>')
            sys.exit()
        elif opt in "-d":
            debris_num = int(arg)
        elif opt in "-s":
            spacecraft_num = int(arg)
        elif opt in "-o":
            output_file = arg

    data["date"] = datetime.datetime.utcnow().isoformat()

    for _ in range(debris_num):
        data["debris"].append(create_debris())

    for _ in range(spacecraft_num):
        data["spacecraft"].append(create_spacecraft())

    if output_file == '':
        with open(ref + "test_input.json", "w") as write_file:
            json.dump(data, write_file)
    else:
        with open(ref + output_file, "w") as write_file:
            json.dump(data, write_file)


if __name__ == '__main__':
    main(sys.argv[1:])
