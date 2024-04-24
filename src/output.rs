use serde::Serialize;

use crate::{bodies, environment};

#[derive(Debug, Serialize)]
pub struct SolarObjectOut {
    pub name: String,    // Name of the solar object
    pub sim_time: f64,   // Simulation time
    pub x_coord: f32,    // Coordinate of object in the x axis
    pub y_coord: f32,    // Coordinate of object in the y axis
    pub z_coord: f32,    // Coordinate of object in the z axis
    pub x_velocity: f64, // Velocity of object in the x axis
    pub y_velocity: f64, // Velocity of object in the y axis
    pub z_velocity: f64, // Velocity of object in the z axis
}

#[derive(Debug, Serialize)]
pub struct SimulationObjectParameters {
    pub id: u32,       // ID of the object
    pub sim_time: f64, // Simulation time
    pub name: String,
    pub soi: String, // Sphere of influence in string form
    pub x_abs_coord: f64,
    /// Absolute coordinate of object in the x axis
    pub y_abs_coord: f64,
    /// Absolute coordinate of object in the y axis
    pub z_abs_coord: f64,
    /// Absolute coordinate of object in the z axis
    pub lat: f64,
    pub long: f64,
    pub altitude: f64,
    pub x_coord: f64,    // Coordinate of object in the x axis
    pub y_coord: f64,    // Coordinate of object in the y axis
    pub z_coord: f64,    // Coordinate of object in the z axis
    pub x_velocity: f64, // Velocity of object in the x axis
    pub y_velocity: f64, // Velocity of object in the y axis
    pub z_velocity: f64, // Velocity of object in the z axis
}

pub trait SimulationOutput {
    fn write_out_object_parameters(&mut self, object_params: SimulationObjectParameters);

    fn write_out_solar_object(&mut self, solar_object: SolarObjectOut);
}

pub mod csv_output {
    use csv;
    use output::{SimulationObjectParameters, SimulationOutput, SolarObjectOut};
    use std::fs;
    use std::path;

    pub struct CSVController {
        object_parameters_writer: csv::Writer<fs::File>,
        solar_object_writer: csv::Writer<fs::File>,
    }

    impl CSVController {
        pub fn new(dir_filepath: &str) -> Self {
            let sub_dirpath = chrono::Utc::now().format("%Y%m%dT%H%M%SZ").to_string();
            let dir_path = path::Path::new(dir_filepath);
            let full_dirpath = dir_path.join(sub_dirpath);
            fs::create_dir_all(full_dirpath.as_path())
                .expect("Cannot create new directory in supplied location.");

            // .unwrap() here is fine as this code is related to initialization
            CSVController {
                object_parameters_writer: csv::Writer::from_path(
                    full_dirpath.join("pose_object_parameters.csv"),
                )
                .unwrap(),
                solar_object_writer: csv::Writer::from_path(
                    full_dirpath.join("pose_solar_objects.csv"),
                )
                .unwrap(),
            }
        }
    }

    impl SimulationOutput for CSVController {
        fn write_out_object_parameters(&mut self, object_params: SimulationObjectParameters) {
            self.object_parameters_writer
                .serialize(object_params)
                .expect(
                    "Failed to write simulation object parameters to the corresponding csv file.",
                );
            // Unwrap here as this is a critical error
            self.object_parameters_writer.flush().unwrap();
        }

        fn write_out_solar_object(&mut self, solar_object: SolarObjectOut) {
            self.solar_object_writer
                .serialize(solar_object)
                .expect("Failed to write simulation solar objects to the corresponding csv file.");
            // Unwrap here as this is a critical error
            self.solar_object_writer.flush().unwrap();
        }
    }
}

pub fn write_out_all_solar_objects(
    env: &environment::Environment,
    output_controller: &mut dyn SimulationOutput,
) {
    output_controller.write_out_solar_object(env.sun_to_output_form());
    output_controller.write_out_solar_object(env.earth_to_output_form());
    output_controller.write_out_solar_object(env.moon_to_output_form());
}

pub fn write_out_all_object_parameters(
    env: &environment::Environment,
    sim_objects: &[bodies::sim_object::SimobjT],
    output_controller: &mut dyn SimulationOutput,
) {
    for sim_obj in sim_objects {
        output_controller.write_out_object_parameters(sim_obj.to_output_form(env.get_sim_time()));
    }
}
