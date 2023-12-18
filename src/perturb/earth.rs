use crate::{
    bodies::{self, sim_object::SimobjT, solar_model::Solarobj},
    environment::Environment,
    output,
    types::{l2_norm, Array3d},
};
use perturb::common;

fn calculate_earth_atmospheric_drag_perturbation(
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    match sim_obj.soi {
        Solarobj::Earth => {}
        _ => return Array3d::default(),
    }

    // return an empty result if sim object earth relative altitude is not within a range which
    // atmospheric drag will be relevant to simulation.
    if !(0.0..=1_000.0 * 1000.0).contains(&sim_obj.coords_fixed.alt) {
        return Array3d::default();
    }

    let atmos_model_result = env
        .earth
        .nrlmsise00_model(env.current_time, &sim_obj.coords_fixed);
    let norm_velocity = l2_norm(&sim_obj.velocity);
    let drag_force = 0.5
        * sim_obj.drag_coeff
        * sim_obj.drag_area
        * atmos_model_result.rho
        * norm_velocity.powf(2.0);
    let drag_force_vector = -drag_force * (sim_obj.velocity / norm_velocity);
    let drag_accel = drag_force_vector / sim_obj.mass;

    if let Some(out) = perturbations_out {
        out.push(output::PerturbationOut {
            id: sim_obj.id,
            sim_time: env.get_sim_time(),
            petrub_type: "drag_earth_nrlmsise00".to_string(),
            acceleration_x_mpss: drag_accel.x,
            acceleration_y_mpss: drag_accel.y,
            acceleration_z_mpss: drag_accel.z,
        })
    }

    drag_accel
}

fn calculate_earth_gravity_perturbation(
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    // Calculate distance between the earth and the sim object using absolute coordinates
    let distance_vector_sim_obj = sim_obj.coords_abs - env.earth.model.state.coords.current_coords;

    // Calculate acceleration due to the earth at the location of the simulation object
    let mut gravity_accel =
        common::newton_gravitational_field(&distance_vector_sim_obj, env.earth.attr.mass);

    // If the simulation object is not in the solar SOI, the gravity field must be compensated
    // This is because perturbation calculations are done relative to the SOI
    gravity_accel = match sim_obj.soi {
        Solarobj::Sun => gravity_accel,
        Solarobj::Earth => {
            bodies::common::ecliptic_to_equatorial(&gravity_accel, env.earth.attr.obliquity)
        }
        Solarobj::Moon => bodies::common::ecliptic_to_equatorial(
            &(gravity_accel
                - common::newton_gravitational_field(
                    &(env.moon.model.state.coords.current_coords
                        - env.earth.model.state.coords.current_coords),
                    env.earth.attr.mass,
                )),
            env.moon.attr.obliquity,
        ),
    };

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
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    calculate_earth_gravity_perturbation(sim_obj, env, perturbations_out)
        + calculate_earth_atmospheric_drag_perturbation(sim_obj, env, perturbations_out)
}
