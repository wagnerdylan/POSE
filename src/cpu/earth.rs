use crate::bodies;
use bodies::KeplerModel;

use super::sim_cpu;

fn calculate_earth_gravity_perturbation(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: &mut Option<&mut Vec<sim_cpu::Perturbation>>,
) -> sim_cpu::PerturbationDelta {
    // Calculate distance between the earth and the sim object using absolute coordinates
    let distance_vector_sim_obj = sim_obj.coords_abs - env.earth.coords.current_coords;

    // Calculate acceleration due to the earth at the location of the simulation object
    let mut gravity_accel = sim_cpu::newton_gravitational_field(
        &distance_vector_sim_obj,
        env.earth.get_solar_object().get_mass_kg(),
    );

    // If the simulation object is not in the solar SOI, the gravity field must be compensated
    // This is because perturbation calculations are done relative to the SOI
    let soi_compensation = match sim_obj.soi {
        bodies::Solarobj::Sun { attr: _ } => None,
        bodies::Solarobj::Earth { attr: _ } => None,
        bodies::Solarobj::Moon { attr: _ } => {
            Some(env.moon.coords.current_coords - env.earth.coords.current_coords)
        }
    };

    // Calculate the difference between the field on simulation object and solar object
    if let Some(soi_info) = soi_compensation {
        gravity_accel = gravity_accel
            - sim_cpu::newton_gravitational_field(
                &soi_info,
                env.earth.get_solar_object().get_mass_kg(),
            );
    }

    let perturbation_delta = sim_cpu::PerturbationDelta {
        id: sim_obj.id,
        sim_time: env.get_sim_time(),
        acceleration: gravity_accel,
    };

    if let Some(out) = perturbations_out {
        out.push(sim_cpu::Perturbation::SolarObject(
            env.earth.get_solar_object().clone(),
            perturbation_delta,
        ))
    }

    perturbation_delta
}

pub fn calculate_earth_perturbations(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: &mut Option<&mut Vec<sim_cpu::Perturbation>>,
) -> sim_cpu::PerturbationDelta {
    // Perturbation due to earth gravity on simulation object
    // TODO combine perturbations if there are more than one
    calculate_earth_gravity_perturbation(sim_obj, env, perturbations_out)
}
