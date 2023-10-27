use serde::Serialize;

use crate::bodies;

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
pub struct PerturbationOut {
    pub id: u32,             // ID of the object perturbation was applied to or calculated for
    pub sim_time: f64,       // Simulation time
    pub petrub_type: String, // Type of perturbation
    pub acceleration_x_mpss: f64, // Acceleration placed on object by perturbing force in x axis
    pub acceleration_y_mpss: f64, // Acceleration placed on object by perturbing force in y axis
    pub acceleration_z_mpss: f64, // Acceleration placed on object by perturbing force in z axis
}

#[derive(Debug, Serialize)]
pub struct SimulationObjectParameters {
    pub id: u32,       // ID of the object
    pub sim_time: f64, // Simulation time
    pub soi: String,   // Sphere of influence in string form
    pub x_abs_coord: f64,
    /// Absolute coordinate of object in the x axis
    pub y_abs_coord: f64,
    /// Absolute coordinate of object in the y axis
    pub z_abs_coord: f64,
    /// Absolute coordinate of object in the z axis
    pub x_coord: f64, // Coordinate of object in the x axis
    pub y_coord: f64,    // Coordinate of object in the y axis
    pub z_coord: f64,    // Coordinate of object in the z axis
    pub x_velocity: f64, // Velocity of object in the x axis
    pub y_velocity: f64, // Velocity of object in the y axis
    pub z_velocity: f64, // Velocity of object in the z axis
}

pub trait SimulationOutput {
    fn write_out_perturbation(&mut self, petrub_out: &PerturbationOut);

    fn write_out_object_parameters(&mut self, object_params: SimulationObjectParameters);

    fn write_out_solar_object(&mut self, solar_object: SolarObjectOut);
}

pub mod csv_output {
    use csv;
    use output::{PerturbationOut, SimulationObjectParameters, SimulationOutput, SolarObjectOut};
    use std::fs;
    use std::path;

    pub struct CSVController {
        perturbation_writer: Option<csv::Writer<fs::File>>,
        object_parameters_writer: csv::Writer<fs::File>,
        solar_object_writer: csv::Writer<fs::File>,
    }

    impl CSVController {
        pub fn new(dir_filepath: &str, should_write_pertub: bool) -> Self {
            let sub_dirpath = chrono::Utc::now().format("%Y%m%dT%H%M%SZ").to_string();
            let dir_path = path::Path::new(dir_filepath);
            let full_dirpath = dir_path.join(sub_dirpath);
            fs::create_dir_all(full_dirpath.as_path())
                .expect("Cannot create new directory in supplied location.");

            // .unwrap() here is fine as this code is related to initialization
            CSVController {
                perturbation_writer: match should_write_pertub {
                    true => Some(
                        csv::Writer::from_path(full_dirpath.join("pose_perturbations.csv"))
                            .unwrap(),
                    ),
                    false => None,
                },
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
        fn write_out_perturbation(&mut self, petrub_out: &PerturbationOut) {
            if let Some(writer) = &mut self.perturbation_writer {
                writer.serialize(petrub_out).expect(
                    "Failed to write simulation perturbation data to the corresponding csv file.",
                );
                // Unwrap here as this is a critical error
                writer.flush().unwrap();
            }
        }

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
    env: &bodies::Environment,
    output_controller: &mut dyn SimulationOutput,
) {
    output_controller.write_out_solar_object(env.sun_to_output_form());
    output_controller.write_out_solar_object(env.earth_to_output_form());
    output_controller.write_out_solar_object(env.moon_to_output_form());
}

pub fn write_out_all_perturbations(
    perturbations: &[PerturbationOut],
    output_controller: &mut dyn SimulationOutput,
) {
    for perturbation in perturbations {
        output_controller.write_out_perturbation(perturbation);
    }
}

pub fn write_out_all_object_parameters(
    env: &bodies::Environment,
    sim_objects: &[bodies::SimobjT],
    output_controller: &mut dyn SimulationOutput,
) {
    for sim_obj in sim_objects {
        output_controller.write_out_object_parameters(sim_obj.to_output_form(env.get_sim_time()));
    }
}
