

pub struct PerturbationOut {
    pub id: u32, // ID of the object perturbation was applied to or calculated for
    pub sim_time: f64, // Simulation time
    pub petrub_type: String, // Type of perturbation
    pub acceleration_x_mpss: f32, // Acceleration placed on object by perturbing force in x axis
    pub acceleration_y_mpss: f32, // Acceleration placed on object by perturbing force in y axis
    pub acceleration_z_mpss: f32, // Acceleration placed on object by perturbing force in z axis
}

pub struct SimulationObjectParameters {
    id: u32, // ID of the object
    sim_time: f64, // Simulation time
    obj_type: String, // Type of Object
    x_coord: f32, // Coordinate of object in the x axis
    y_coord: f32, // Coordinate of object in the y axis
    z_coord: f32, // Coordinate of object in the z axis
    x_velocity: f32, // Velocity of object in the x axis
    y_velocity: f32, // Velocity of object in the y axis
    z_velocity: f32, // Velocity of object in the z axis
}

pub trait SimulationOutput {

    fn write_out_perturbation(&self, petrub_out: PerturbationOut);

    fn write_out_object_parameters(&self, object_params: SimulationObjectParameters);

}

pub mod csv_output {
    use output::{SimulationOutput, PerturbationOut, SimulationObjectParameters};

    pub struct CSVController {
        // TODO add csv specific stuff
    }

    impl CSVController {

        pub fn new(dir_path: &str) -> Self {
            CSVController{}
        }

    }

    impl SimulationOutput for CSVController{
        fn write_out_perturbation(&self, petrub_out: PerturbationOut) {
            unimplemented!()
        }

        fn write_out_object_parameters(&self, object_params: SimulationObjectParameters) {
            unimplemented!()
        }
    }
}