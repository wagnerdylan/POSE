use crate::bodies;
use bodies::KeplerModel;

use super::sim_cpu;

fn calculate_moon_gravity_perturbation(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: &mut Option<&mut Vec<sim_cpu::Perturbation>>,
) -> sim_cpu::PerturbationDelta {
    // Calculate distance between the moon and the sim object using absolute coordinates
    let distance_vector_sim_obj = sim_obj.coords_abs - env.current_moon_coords;

    // Calculate acceleration due to the moon at the location of the simulation object
    let gravity_accel = sim_cpu::newton_gravitational_field(
        &distance_vector_sim_obj,
        env.moon.get_solar_object().get_mass_kg(),
    );

    let perturbation_delta = sim_cpu::PerturbationDelta {
        id: sim_obj.id,
        sim_time: env.get_sim_time(),
        acceleration: gravity_accel,
    };

    if let Some(out) = perturbations_out {
        out.push(sim_cpu::Perturbation::SolarObject(
            env.moon.get_solar_object().clone(),
            perturbation_delta,
        ))
    }

    perturbation_delta
}

pub fn calculate_moon_perturbations(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: &mut Option<&mut Vec<sim_cpu::Perturbation>>,
) -> sim_cpu::PerturbationDelta {
    // Perturbation due to moon gravity on simulation object
    let gravity_perturbation = calculate_moon_gravity_perturbation(sim_obj, env, perturbations_out);

    // TODO combine perturbations if there are more than one
    gravity_perturbation
}
