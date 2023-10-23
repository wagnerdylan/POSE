use super::bodies;

use chrono::DateTime;
use clap::Parser;
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
        long,
        help = "value to be used as the simulation time step in seconds.",
        default_value_t = 1.0
    )]
    pub sim_time_step: f32,
    #[arg(long, help = "value to be used as the update period for solar objects.", default_value_t = 3600.0 * 12.0 )]
    pub sim_solar_step: f32,
}

/// Main entry point into the init sequence
///
/// ### Argument
/// * 'file' - The name of the input file containing the bodies
///
/// ### Return
///      A vector of bodies from the input file.
///      The datetime delta from year 2000-01-01
///
pub fn parse_input(file: &str) -> (Vec<bodies::SimobjT>, DateTime<chrono::Utc>) {
    let mut sim_bodies: Vec<bodies::SimobjT> = Vec::new();

    let ser_objs = read_object_from_file(file).unwrap();

    //add objects to sim_bodies
    for mut elem in ser_objs.debris {
        elem.sim_object_type = bodies::SimObjectType::Debris;
        sim_bodies.push(elem);
    }

    for mut elem in ser_objs.spacecraft {
        elem.sim_object_type = bodies::SimObjectType::Spacecraft;
        sim_bodies.push(elem);
    }

    assign_id(&mut sim_bodies);

    let datetime = ser_objs.date;
    let datetime_obj = datetime
        .parse::<DateTime<chrono::Utc>>()
        .expect("Input file contains invalid datetime format, expected ISO 8601 format.");

    (sim_bodies, datetime_obj)
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
fn read_object_from_file<P: AsRef<Path>>(path: P) -> Result<bodies::InitData, Box<dyn Error>> {
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
fn assign_id(sim_bodies: &mut Vec<bodies::SimobjT>) {
    let mut id_inc: u32 = 1;

    for body in sim_bodies {
        body.id = id_inc;
        id_inc += 1;
    }
}
