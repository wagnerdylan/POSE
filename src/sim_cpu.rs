use crate::bodies;
use crate::output;
use input::SimulationParameters;
use std::string::ToString;
use strum_macros::Display;

pub struct PerturbationDelta {
    id: u32,
    sim_time: f64,
    acceleration_x_mpss: f32,
    acceleration_y_mpss: f32,
    acceleration_z_mpss: f32,
}

impl Default for PerturbationDelta {
    fn default() -> Self {
        Self {
            id: 0,
            sim_time: 0.0,
            acceleration_x_mpss: 0.0,
            acceleration_y_mpss: 0.0,
            acceleration_z_mpss: 0.0,
        }
    }
}

#[derive(Display)]
pub enum Perturbation {
    #[strum(serialize = "solar_obj")]
    SolarObject(bodies::Solarobj, PerturbationDelta),
}

impl Perturbation {
    fn into_output_form(self) -> output::PerturbationOut {
        output::PerturbationOut {
            id: match &self {
                Perturbation::SolarObject(_, perturb_delta) => perturb_delta.id,
            },
            sim_time: match &self {
                Perturbation::SolarObject(_, perturb_delta) => perturb_delta.sim_time,
            },
            petrub_type: match &self {
                Perturbation::SolarObject(solar_obj, _) => {
                    format!("{}_{}", self.to_string(), solar_obj.to_string())
                }
            },
            acceleration_x_mpss: match &self {
                Perturbation::SolarObject(_, perturb_delta) => perturb_delta.acceleration_x_mpss,
            },
            acceleration_y_mpss: match &self {
                Perturbation::SolarObject(_, perturb_delta) => perturb_delta.acceleration_y_mpss,
            },
            acceleration_z_mpss: match &self {
                Perturbation::SolarObject(_, perturb_delta) => perturb_delta.acceleration_z_mpss,
            },
        }
    }
}

/// Write out all perturbations to output controller.
///
/// ### Arguments
/// * 'perturbations' - Vector containing perturbations to be written out.
/// * 'output_controller' - Controller object used to facilitate the output of perturbation data.
///
fn write_out_all_perturbations(
    perturbations: Vec<Perturbation>,
    output_controller: &mut Box<dyn output::SimulationOutput>,
) {
    for perturbation in perturbations {
        output_controller.write_out_perturbation(perturbation.into_output_form());
    }
}

fn write_out_all_object_parameters(
    _sim_objects: &Vec<bodies::SimobjT>,
    _output_controller: &mut Box<dyn output::SimulationOutput>,
) {
    unimplemented!();
}

fn write_out_all_solar_objects(
    env: &bodies::Environment,
    output_controller: &mut Box<dyn output::SimulationOutput>,
) {
    for solar_object in env.get_solar_objects() {
        output_controller.write_out_solar_object(solar_object.to_output_form(env.sim_time_s));
    }
}

fn l2_norm(x: ndarray::ArrayView1<f64>) -> f64 {
    x.dot(&x).sqrt()
}

fn normalize(mut x: ndarray::Array1<f64>) -> ndarray::Array1<f64> {
    let norm = l2_norm(x.view());
    x.mapv_inplace(|e| e/norm);
    x
}

/// Module used to apply perturbation calculations on individual bodies
mod cowell_perturb {
    use crate::bodies;
    use crate::sim_cpu::{Perturbation, PerturbationDelta};

    /// Apply all perturbations handled by POSE. This includes:
    /// * 'Solar Body Earth'
    /// * 'Solar Body Moon'
    /// * 'Solar Body Sun'
    /// TODO add more
    ///
    /// ### Parameters
    /// * 'sim_obj' - The object basis for calculation and apply
    /// * 'env' - The Simulation environment
    /// * 'do_return_peturb' - true if vector should be returned, false otherwise
    ///
    /// ### Return
    ///     A vector of perturbation deltas in _do_return_peturb is true or none.
    ///
    pub fn apply_perturbations(
        _sim_obj: &mut bodies::SimobjT,
        _env: &bodies::Environment,
        _do_return_peturb: bool,
    ) -> Option<Vec<Perturbation>> {
        // TODO collect a vector of all perturbations
        // TODO sum perturbations
        // TODO apply sum of perturbations to sim_obj
        // TODO emit vector of all perturbations

        unimplemented!();
    }


    /// Calculate perturbations due to solar system objects.
    ///
    /// ### Parameters
    /// * 'sim_obj' - The object basis for calculation
    /// * 'env' - The Simulation environment
    /// * 'do_return_peturb' - true if vector should be returned, false otherwise
    ///
    /// ### Return
    ///     A struct of size two containing
    ///         (Total perturbation delta, individual perturbation deltas or none)
    ///
    fn calc_planet_perturb(
        sim_obj: &bodies::SimobjT,
        env: &bodies::Environment,
        _do_return_peturb: bool,
    ) -> (PerturbationDelta, Option<Vec<Perturbation>>) {

        // TODO function which uses vector form of Newton's law of universal gravitation
        fn newton_gravitation(sim_obj: &bodies::SimobjT, env: &bodies::Environment)
            -> ndarray::Array1<f64>
        {
            unimplemented!();
            // TODO Calculate L2 Norm from sim_obj to env centric
            // TODO Calculate unit vector for perturbation
            // TODO Calculate force using Newton's law of universal gravitation
            // TODO Return the force vector
        }

        unimplemented!();
    }
}

/// Main entry point into the cpu_sim module, gathers all needed data for orbit modeling
/// using Cowell's method.
pub fn simulate(
    _sim_bodies: Vec<bodies::SimobjT>,
    mut env: bodies::Environment,
    mut output_controller: Box<dyn output::SimulationOutput>,
    sim_params: SimulationParameters,
) {
    loop {

        // Update solar objs
        if env.sim_time_s > env.last_day_update_s + sim_params.sim_solar_step as f64 {
            write_out_all_solar_objects(&env, &mut output_controller);
            env.update();
        }

        // Move forward simulation by step
        env.sim_time_s += sim_params.sim_time_step as f64;
    }
}
