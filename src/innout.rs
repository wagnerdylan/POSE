
extern crate serde_json;
extern crate chrono;

use super::bodies;

use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use self::chrono::DateTime;
use serde::de::value::StrDeserializer;
use self::chrono::offset::TimeZone;

/// Main entry point into the init sequence
///
/// ### Argument
/// * 'file' - The name of the input file containing the bodies
///
/// ### Return
///      A vector of bodies from the input file.
///      The datetime delta from year 2000-01-01
///
pub fn parse_inpt(file: &str) -> (Vec<bodies::SimobjT>, f32){
    let mut sim_bodies: Vec<bodies::SimobjT> = Vec::new();

    let ser_objs = read_object_from_file(file).unwrap();

    //add objects to sim_bodies
    for elem in ser_objs.debris {
        //println!("id {}", elem.type_of());
        let p = Box::new(elem.clone());
        sim_bodies.push(p);
    }

    for elem in ser_objs.spacecraft {
        //println!("id {}", elem.type_of());
        let p = Box::new(elem.clone());
        sim_bodies.push(p);
    }

    assign_id(&mut sim_bodies);

    let datetime = ser_objs.date;
    let dt_delta = datetime_to_days(datetime);

    (sim_bodies, dt_delta)
}

/// Calculates a delta for provided datetime from 0/Jan/2000 00:00 UTC
///
/// ### Argument
/// * 'datetime' - User provided datetime string.
///
/// ### Return
///     The delta from 0/Jan/2000 00:00 UTC in days.
///
fn datetime_to_days(datetime: String) -> f32 {

    let origin_dt = chrono::Utc.ymd(2000,1,1).and_hms(0,0,0);
    let datetime_obj= datetime.parse::<DateTime<chrono::Utc>>().expect(
        "Input file contains invalid datetime format, expected ISO 8601 format.");

    1.15741e-5f32 * (datetime_obj- origin_dt).num_seconds() as f32

}

/// Function responsable for handling opeing the file and conneting the
/// serde reader.
///
/// ### Argument
/// * 'path' - The path to the input bodies json file.
///
/// ### Return
///      A result object loaded with an IO error on failure or the serde reader on
///      success.
///
fn read_object_from_file<P: AsRef<Path>>(path: P) -> Result<bodies::InitData, Box<Error>> {
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
        *body.id_mut() = id_inc;
        id_inc += 1;
    }
}