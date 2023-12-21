use crate::{
    bodies::{common::ecliptic_to_equatorial, sim_object::SimobjT, solar_model::Solarobj},
    environment::Environment,
    output,
    types::{l2_norm, Array3d},
};

use super::common::{newton_gravitational_field, newton_gravitational_field_third_body};

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
    let (r_0, r_tb, ob_0) = match sim_obj.soi {
        Solarobj::Sun => (
            env.sun.model.state.coords.current_coords,
            Some(env.earth.model.state.coords.current_coords),
            env.sun.attr.obliquity,
        ),
        Solarobj::Earth => (
            env.earth.model.state.coords.current_coords,
            None,
            env.earth.attr.obliquity,
        ),
        Solarobj::Moon => (
            env.moon.model.state.coords.current_coords,
            Some(env.earth.model.state.coords.current_coords),
            env.moon.attr.obliquity,
        ),
    };

    let accel_ecliptic = if r_tb.is_some() {
        newton_gravitational_field_third_body(
            &sim_obj.coords_abs,
            &r_0,
            &r_tb.unwrap(),
            env.earth.attr.mass,
        )
    } else {
        newton_gravitational_field(&sim_obj.coords_abs, &r_0, env.earth.attr.mass)
    };

    let accel = ecliptic_to_equatorial(&accel_ecliptic, ob_0);

    if let Some(out) = perturbations_out {
        out.push(output::PerturbationOut {
            id: sim_obj.id,
            sim_time: env.get_sim_time(),
            petrub_type: "solar_obj_earth".to_string(),
            acceleration_x_mpss: accel.x,
            acceleration_y_mpss: accel.y,
            acceleration_z_mpss: accel.z,
        })
    }

    accel
}

pub fn calculate_earth_perturbations(
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    calculate_earth_gravity_perturbation(sim_obj, env, perturbations_out)
        + calculate_earth_atmospheric_drag_perturbation(sim_obj, env, perturbations_out)
}
