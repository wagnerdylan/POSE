use bodies::sim_object::{SimObjectType, SimobjT};
use chrono::{DateTime, Utc};
use clap::Parser;
use serde::Deserialize;
use std::{error::Error, fs::File, io::BufReader, path::Path};

use crate::bodies::sim_object::PerturbationStore;

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
        default_value_t = 0.1
    )]
    pub write_period: f64,
    #[arg(
        long,
        help = "value to be used as the simulation time step in seconds.",
        default_value_t = 1.0
    )]
    pub sim_time_step: f32,
    #[arg(
        long,
        help = "value to be used as the update period for solar objects.",
        default_value_t = 60.0 * 60.0
    )]
    pub sim_solar_step: f32,
    #[arg(
        long,
        help = "period which to run collision detection.",
        default_value_t = 1.0 * 60.0
    )]
    pub collision_check_period: f32,
    #[arg(
        long,
        help = "flag used to indicate perturbations should be written out.",
        default_value_t = false
    )]
    pub write_perturbations: bool,
}

#[derive(Deserialize, Clone, Default, Debug)]
pub struct SwIndex(
    pub DateTime<chrono::Utc>, // Index start time.
    pub DateTime<chrono::Utc>, // Index end time.
    pub f64,                   // AP value.
    pub f64,                   // F10.7_OBS value.
    pub f64,                   // F10.7_OBS_CENTER8 value.
);

impl SwIndex {
    #[inline]
    pub(crate) fn cmp(&self, query_time: DateTime<Utc>) -> std::cmp::Ordering {
        if self.0.le(&query_time) && query_time.le(&self.1) {
            std::cmp::Ordering::Equal
        } else if query_time.lt(&self.0) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Less
        }
    }
}

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
    pub sim_time_step: f32,
    pub sim_solar_step: f32,
    pub collision_check_period: f32,
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

    if sim_params.write_perturbations {
        init_perturbation_objects(&mut sim_bodies);
    }

    let runtime_params = RuntimeParameters {
        date: ser_objs.date,
        halt_date: ser_objs.halt_date,
        write_period: sim_params.write_period,
        sim_time_step: sim_params.sim_time_step,
        sim_solar_step: sim_params.sim_solar_step,
        collision_check_period: sim_params.collision_check_period,
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

/// Initialize perturbation store object within simulation body to
/// Some() value. This indicates simulation output should include per
/// body perturbation accelerations.
///
/// ### Argument
/// * 'sim_bodies' - A vector containing both debris and spacecraft objects.
///
fn init_perturbation_objects(sim_bodies: &mut Vec<SimobjT>) {
    for body in sim_bodies {
        body.perturb_store = Some(PerturbationStore::default())
    }
}
