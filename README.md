![POSE Logo](/images/POSE_Logo.svg)\
Parallel Orbital Simulation Environment

# Usage
pose [OPTIONS] \<INPUT\>

# Building
llvm must be installed on the build system in order to build POSE. On Windows, llvm may be installed via choco using the following command: 
`choco install llvm`

On install of llvm (if needed) pose may be built using cargo such as the following: `cargo build --release`

# Basic Sanity Checking
Under data/, iss.json is used for basic sanity checks on simulation correctness. ISS trajectory data is gathered from https://spotthestation.nasa.gov/trajectory_data.cfm provided by NASA.