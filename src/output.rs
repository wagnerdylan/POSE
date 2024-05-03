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
    /// Heliocentric coordinate of object in the x axis
    pub x_coord_helio: f64,
    /// Heliocentric coordinate of object in the y axis
    pub y_coord_helio: f64,
    /// Heliocentric coordinate of object in the z axis
    pub z_coord_helio: f64,
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

#[derive(Debug, Serialize)]
pub struct SimulationObjectPerturbationOut {
    pub id: u32,       // ID of the object
    pub sim_time: f64, // Simulation time
    pub name: String,
    pub perturb_type: String,
    pub x_accel: f64,
    pub y_accel: f64,
    pub z_accel: f64,
}

pub trait SimulationOutput {
    fn write_out_object_parameters(&mut self, object_params: SimulationObjectParameters);

    fn write_out_object_perturbation(&mut self, perturb_object: SimulationObjectPerturbationOut);

    fn write_out_solar_object(&mut self, solar_object: SolarObjectOut);
}

pub mod csv_output {
    use csv;
    use output::{SimulationObjectParameters, SimulationOutput, SolarObjectOut};
    use std::fs;
    use std::path;

    pub struct CSVController {
        object_parameters_writer: csv::Writer<fs::File>,
        object_perturbation_writer: csv::Writer<fs::File>,
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
                object_perturbation_writer: csv::Writer::from_path(
                    full_dirpath.join("pose_object_perturbations.csv"),
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
            // Treat as critical error.
            self.object_parameters_writer.flush().unwrap();
        }

        fn write_out_solar_object(&mut self, solar_object: SolarObjectOut) {
            self.solar_object_writer
                .serialize(solar_object)
                .expect("Failed to write simulation solar objects to the corresponding csv file.");
            // Treat as critical error.
            self.solar_object_writer.flush().unwrap();
        }

        fn write_out_object_perturbation(
            &mut self,
            perturb_object: super::SimulationObjectPerturbationOut,
        ) {
            self.object_perturbation_writer
                .serialize(perturb_object)
                .expect(
                "Failed to write out simulation object perturbation to the corresponding csv file.",
            );
            // Treat as critical error.
            self.object_perturbation_writer.flush().unwrap();
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

pub fn write_out_all_object_perturbations(
    env: &environment::Environment,
    sim_objects: &[bodies::sim_object::SimobjT],
    output_controller: &mut dyn SimulationOutput,
) {
    for sim_obj in sim_objects {
        if let Some(store) = &sim_obj.perturb_store {
            for perturb_out in store.to_output_form(sim_obj.id, &sim_obj.name, env.get_sim_time()) {
                output_controller.write_out_object_perturbation(perturb_out);
            }
        }
    }
}
