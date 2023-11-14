use crate::{bodies, output, types::Array3d};
use bodies::KeplerModel;
use perturb::common;

fn calculate_moon_gravity_perturbation(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    // Calculate distance between the moon and the sim object using absolute coordinates
    let distance_vector_sim_obj = sim_obj.coords_abs - env.moon.coords.current_coords;

    // Calculate acceleration due to the moon at the location of the simulation object
    let gravity_accel = common::newton_gravitational_field(
        &distance_vector_sim_obj,
        env.moon.get_solar_object().get_mass_kg(),
    );

    if let Some(out) = perturbations_out {
        out.push(output::PerturbationOut {
            id: sim_obj.id,
            sim_time: env.get_sim_time(),
            petrub_type: "solar_obj_moon".to_string(),
            acceleration_x_mpss: gravity_accel.x,
            acceleration_y_mpss: gravity_accel.y,
            acceleration_z_mpss: gravity_accel.z,
        })
    }

    gravity_accel
}

pub fn calculate_moon_perturbations(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    calculate_moon_gravity_perturbation(sim_obj, env, perturbations_out)
}
