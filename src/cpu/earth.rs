use super::sim_cpu::{self};
use crate::{bodies, output, types::Array3d};
use bodies::KeplerModel;

fn calculate_earth_atmospheric_drag_perturbation(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    // TODO stub function

    if let Some(out) = perturbations_out {
        out.push(output::PerturbationOut {
            id: sim_obj.id,
            sim_time: env.get_sim_time(),
            petrub_type: "drag_earth".to_string(),
            acceleration_x_mpss: 0.0,
            acceleration_y_mpss: 0.0,
            acceleration_z_mpss: 0.0,
        })
    }

    Array3d::default()
}

fn calculate_earth_gravity_perturbation(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
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

    if let Some(out) = perturbations_out {
        out.push(output::PerturbationOut {
            id: sim_obj.id,
            sim_time: env.get_sim_time(),
            petrub_type: "solar_obj_earth".to_string(),
            acceleration_x_mpss: gravity_accel.x,
            acceleration_y_mpss: gravity_accel.y,
            acceleration_z_mpss: gravity_accel.z,
        })
    }

    gravity_accel
}

pub fn calculate_earth_perturbations(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    calculate_earth_gravity_perturbation(sim_obj, env, perturbations_out)
        + calculate_earth_atmospheric_drag_perturbation(sim_obj, env, perturbations_out)
}
