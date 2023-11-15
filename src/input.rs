use bodies::sim_object::{SimObjectType, SimobjT};
use chrono::DateTime;
use clap::Parser;
use serde::Deserialize;
use std::{error::Error, fs::File, io::BufReader, path::Path};

#[derive(Parser, Debug)]
#[command(
    version,
    name = "Parallel Orbital Simulation Environment (POSE)",
    long_about = "Simulation aimed to model the orbital environment around Earth for bodies at all magnitudes."
)]
pub struct SimulationParameters {
    #[arg(
        short,
        long,
        help = "json file containing information on bodies at initialization."
    )]
    pub configuration: String,
    #[arg(short, long, help = "name of desired output directory.")]
    pub output_dir: String,
    #[arg(
        short,
        long,
        help = "simulation duration between output writes. This parameter should be expressed in seconds.",
        default_value_t = 1.0
    )]
    pub write_period: f64,
    #[arg(
        long,
        help = "flag used to specify perturbation results per simulation objects.",
        default_value_t = false
    )]
    pub write_out_pertub: bool,
    #[arg(
        long,
        help = "value to be used as the simulation time step in seconds.",
        default_value_t = 1.0
    )]
    pub sim_time_step: f32,
    #[arg(
        long,
        help = "value to be used as the update period for solar objects.",
        default_value_t = 3600.0
    )]
    pub sim_solar_step: f32,
}

#[derive(Deserialize)]
pub struct SwIndex(DateTime<chrono::Utc>, DateTime<chrono::Utc>, f64, f64, f64);

#[derive(Deserialize)]
pub struct InitData {
    pub date: DateTime<chrono::Utc>,      // Datetime in ISO 8601 format.
    pub halt_date: DateTime<chrono::Utc>, // Date used as simulation stopping condition.
    pub spacecraft: Vec<SimobjT>,         // Spacecraft objects.
    pub debris: Vec<SimobjT>,             // Debris objects.
    pub earth_sw: Vec<SwIndex>,           // Space Weather Data for Earth atmospheric model.
}

pub struct EnvInitData {
    pub earth_sw: Vec<SwIndex>, // Space Weather Data for Earth atmospheric model.
}

#[derive(Default)]
pub struct RuntimeParameters {
    pub date: DateTime<chrono::Utc>,      // Datetime in ISO 8601 format.
    pub halt_date: DateTime<chrono::Utc>, // Date used as simulation stopping condition.
    pub write_period: f64,
    pub write_out_pertub: bool,
    pub sim_time_step: f32,
    pub sim_solar_step: f32,
}

pub fn collect_simulation_inputs(
    sim_params: &SimulationParameters,
) -> (Vec<SimobjT>, RuntimeParameters, EnvInitData) {
    let mut sim_bodies: Vec<SimobjT> = Vec::new();

    let ser_objs = read_object_from_file(&sim_params.configuration).unwrap();

    //add objects to sim_bodies
    for mut elem in ser_objs.debris {
        elem.sim_object_type = SimObjectType::Debris;
        sim_bodies.push(elem);
    }

    for mut elem in ser_objs.spacecraft {
        elem.sim_object_type = SimObjectType::Spacecraft;
        sim_bodies.push(elem);
    }

    assign_id(&mut sim_bodies);

    let runtime_params = RuntimeParameters {
        date: ser_objs.date,
        halt_date: ser_objs.halt_date,
        write_period: sim_params.write_period,
        write_out_pertub: sim_params.write_out_pertub,
        sim_time_step: sim_params.sim_time_step,
        sim_solar_step: sim_params.sim_solar_step,
    };

    let env_init_data = EnvInitData {
        earth_sw: ser_objs.earth_sw,
    };

    (sim_bodies, runtime_params, env_init_data)
}

/// Function responsible for handling opening the file and connecting the
/// serde reader.
///
/// ### Argument
/// * 'path' - The path to the input bodies json file.
///
/// ### Return
///      A result object loaded with an IO error on failure or the serde reader on
///      success.
///
fn read_object_from_file<P: AsRef<Path>>(path: &P) -> Result<InitData, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `Objects`.
    let u = serde_json::from_reader(reader)?;

    // Return the `Objects`.
    Ok(u)
}

/// Adds an sequential id value to each of the simulation bodies.
///
/// ### Argument
/// * 'sim_bodies' - A vector containing both debris and spacecraft objects.
///
fn assign_id(sim_bodies: &mut Vec<SimobjT>) {
    let mut id_inc: u32 = 1;

    for body in sim_bodies {
        body.id = id_inc;
        id_inc += 1;
    }
}
