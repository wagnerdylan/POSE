use crate::output;
use crate::{bodies, input};
use std::string::ToString;
use strum_macros::Display;
use types::Array3d;

use super::sol;
use super::{earth, moon};

// Gravitational constant 6.674×10−11
const G: f64 = 6.674e-11;
const MAX_NUM_OF_PERTURBATIONS: usize = 3;

#[derive(Clone, Copy)]
pub struct PerturbationDelta {
    pub(crate) id: u32,
    pub(crate) sim_time: f64,
    pub(crate) acceleration: Array3d,
}

impl Default for PerturbationDelta {
    fn default() -> Self {
        Self {
            id: 0,
            sim_time: 0.0,
            acceleration: Array3d::default(),
        }
    }
}

#[derive(Display)]
pub enum Perturbation {
    #[strum(serialize = "solar_obj")]
    SolarObject(bodies::Solarobj, PerturbationDelta),
}

impl Perturbation {
    pub fn to_output_form(&self) -> output::PerturbationOut {
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
                Perturbation::SolarObject(_, perturb_delta) => perturb_delta.acceleration.x,
            },
            acceleration_y_mpss: match &self {
                Perturbation::SolarObject(_, perturb_delta) => perturb_delta.acceleration.y,
            },
            acceleration_z_mpss: match &self {
                Perturbation::SolarObject(_, perturb_delta) => perturb_delta.acceleration.z,
            },
        }
    }
}

pub fn l2_norm(x: &Array3d) -> f64 {
    x.dot(x).sqrt()
}

pub fn normalize(x: &Array3d, l2_norm_precalc: Option<f64>) -> Array3d {
    let norm = match l2_norm_precalc {
        Some(val) => val,
        None => l2_norm(x),
    };

    Array3d {
        x: x.x / norm,
        y: x.y / norm,
        z: x.z / norm,
    }
}

pub fn newton_gravitational_field(distance_vector: &Array3d, planet_mass_kg: f64) -> Array3d {
    let l2_dist = l2_norm(distance_vector);
    // Calculate unit vector for perturbation
    let unit_vector = normalize(distance_vector, Some(l2_dist));
    // Calculate force using Newton's law of universal gravitation
    unit_vector * (-G * (planet_mass_kg / l2_dist.powi(2)))
}

pub fn apply_perturbations(
    sim_obj: &mut bodies::SimobjT,
    env: &bodies::Environment,
    step_time_s: f64,
    perturbations_out: &mut Option<&mut Vec<Perturbation>>,
) {
    sim_obj.coords_abs = env.calculate_abs_coords(sim_obj);

    let mut net_acceleration = Array3d::default();
    // Calculate the pertubation forces for all planetary objects
    net_acceleration = net_acceleration
        + sol::calculate_solar_perturbations(sim_obj, env, perturbations_out).acceleration;
    net_acceleration = net_acceleration
        + earth::calculate_earth_perturbations(sim_obj, env, perturbations_out).acceleration;
    net_acceleration = net_acceleration
        + moon::calculate_moon_perturbations(sim_obj, env, perturbations_out).acceleration;

    // Calculate the velocity change from net acceleration
    let velocity_delta = net_acceleration * step_time_s;
    // Calculate new velocity for the given simulation object
    let updated_sim_obj_velocity = velocity_delta + sim_obj.velocity;

    // Calculate the position change from the updated velocity
    let position_delta = updated_sim_obj_velocity * step_time_s;
    // Calculate the new position for the simulation object
    let updated_sim_obj_coords = position_delta + sim_obj.coords;

    // Update the new values within the simulation object
    sim_obj.velocity = updated_sim_obj_velocity;
    sim_obj.coords = updated_sim_obj_coords;
}

/// Main entry point into the cpu_sim module, gathers all needed data for orbit modeling
/// using Cowell's method.
pub fn simulate(
    mut sim_bodies: Vec<bodies::SimobjT>,
    mut env: bodies::Environment,
    mut output_controller: Box<dyn output::SimulationOutput>,
    sim_params: input::SimulationParameters,
) {
    // Allocate large buffer for holding perturbations
    let mut perturbation_vec: Vec<Perturbation> =
        Vec::with_capacity(sim_bodies.len() * MAX_NUM_OF_PERTURBATIONS);

    loop {
        // Calculate and apply perturbations for every object
        // TODO parallelize this
        for sim_obj in sim_bodies.iter_mut() {
            apply_perturbations(
                sim_obj,
                &env,
                sim_params.sim_time_step as f64,
                &mut Some(&mut perturbation_vec),
            );
        }

        output::write_out_all_object_parameters(&env, &sim_bodies, output_controller.as_mut());
        output::write_out_all_solar_objects(&env, output_controller.as_mut());
        output::write_out_all_perturbations(&mut perturbation_vec, output_controller.as_mut());

        env.advance_simulation_environment(&sim_params);

        perturbation_vec.clear();
    }
}
