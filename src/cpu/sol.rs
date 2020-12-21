use crate::bodies;
use bodies::KeplerModel;

use super::sim_cpu;

fn calculate_solar_gravity_perturbation(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: Option<&mut Vec<sim_cpu::Perturbation>>,
) -> sim_cpu::PerturbationDelta {
    // Calculate distance between the sun and the sim object using absolute coordinates
    let distance_vector_sim_obj = env.sun.distance_between_sim_object(sim_obj);

    // Calculate acceleration due to the sun at the location of the simulation object
    let mut gravity_accel = sim_cpu::newton_gravitational_field(
        &distance_vector_sim_obj,
        env.sun.get_solar_object().get_mass_kg(),
    );

    // If the simulation object is not in the solar SOI, the gravity field must be compensated
    // This is because perturbation calculations are done relative to the SOI
    let soi_compensation = match sim_obj.soi {
        bodies::Solarobj::Sun { attr: _ } => None,
        bodies::Solarobj::Earth { attr: _ } => Some((
            env.earth.get_solar_object().get_mass_kg(),
            env.current_earth_coords,
        )),
        bodies::Solarobj::Moon { attr: _ } => Some((
            env.moon.get_solar_object().get_mass_kg(),
            env.current_moon_coords,
        )),
    };

    // Calculate the difference between the field on simulation object and solar object
    if let Some(soi_info) = soi_compensation {
        gravity_accel =
            gravity_accel - sim_cpu::newton_gravitational_field(&soi_info.1, soi_info.0);
    }

    let perturbation_delta = sim_cpu::PerturbationDelta {
        id: sim_obj.id,
        sim_time: env.get_sim_time(),
        acceleration: gravity_accel,
    };

    if let Some(out) = perturbations_out {
        out.push(sim_cpu::Perturbation::SolarObject(
            env.sun.get_solar_object().clone(),
            perturbation_delta,
        ))
    }

    perturbation_delta
}

/// Calculate all solar perturbations on sim_obj
///
/// ### Arguments
/// * 'sim_obj' - Simulation object
/// * 'env' - Simulation environment
/// * 'perturbations_out' - Array of individual perturbations wrapped in a option.
///    This array should be appended to if option value is some
/// * 'step_time_s' - Simulation step time
///
/// ### Returns
/// Combined perturbations
///
pub fn calculate_solar_perturbations(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: Option<&mut Vec<sim_cpu::Perturbation>>,
) -> sim_cpu::PerturbationDelta {
    // Perturbation due to solar gravity on simulation object
    let gravity_perturbation =
        calculate_solar_gravity_perturbation(sim_obj, env, perturbations_out);

    // TODO combine perturbations if there are more than one
    gravity_perturbation
}
