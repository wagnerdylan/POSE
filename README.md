![POSE Logo](/images/POSE_Logo.svg)\
Parallel Orbital Simulation Environment

## Motivation
The main goal of POSE is to provide a means to understand how space debris impacts operation of spacecraft across lengthy periods of time. Ancillary goals include education and inspiration of orbit simulation to the general public through public facing interfaces into POSE.

## Configuration and Usage
POSE is configured via a single JSON configuration file which is crafted via hand and generation programs provided under this repo in `util`. 

TODO add more info.

## Building
llvm must be installed on the build system in order to build POSE. On Windows, llvm may be installed via choco using the following command: 
`choco install llvm`

NASA C-SPICE must also be installed on the host machine following the instructions provided [here](https://naif.jpl.nasa.gov/naif/toolkit_C.html). The environment variable
`CSPICE_DIR` must be set to the top-level directory of the SPICE directory in-order to pull in required SPICE libs.

On install of llvm (if needed) pose may be built using cargo such as the following: `cargo build --release`

## Basic Sanity Checking
Under data/, iss.json is used for basic sanity checks on simulation correctness. ISS trajectory data is gathered from https://spotthestation.nasa.gov/trajectory_data.cfm provided by NASA.