use crate::{
    bodies::{
        common::ecliptic_to_equatorial,
        sim_object::{PerturbationDefinition, SimobjT},
        solar_model::Solarobj,
    },
    environment::Environment,
    types::{l2_norm, Array3d},
};

use super::common::{newton_gravitational_field, newton_gravitational_field_third_body};

fn calculate_earth_atmospheric_drag_perturbation(
    sim_obj: &mut SimobjT,
    env: &Environment,
) -> Array3d {
    match sim_obj.state.soi {
        Solarobj::Earth => {}
        _ => return Array3d::default(),
    }

    // return an empty result if sim object earth relative altitude is not within a range which
    // atmospheric drag will be relevant to simulation.
    if !(0.0..=1_000.0 * 1000.0).contains(&sim_obj.state.coords_fixed.alt) {
        if let Some(perturb_store) = &mut sim_obj.perturb_store {
            perturb_store.earth_drag = None;
        }
        return Array3d::default();
    }

    let atmos_model_result = env
        .earth
        .nrlmsise00_model(env.current_time, &sim_obj.state.coords_fixed);
    let norm_velocity = l2_norm(&sim_obj.state.velocity);
    let drag_force = 0.5
        * sim_obj.drag_coeff
        * sim_obj.drag_area
        * atmos_model_result.rho
        * norm_velocity.powf(2.0);
    let drag_force_vector = -drag_force * (sim_obj.state.velocity / norm_velocity);

    let perturb = drag_force_vector / sim_obj.mass;

    if let Some(perturb_store) = &mut sim_obj.perturb_store {
        perturb_store.earth_drag = Some(PerturbationDefinition {
            perturb_name: "earth_drag_nrlmsise00".to_string(),
            x_accel: perturb.x,
            y_accel: perturb.y,
            z_accel: perturb.z,
        });
    }

    perturb
}

fn calculate_earth_gravity_perturbation(sim_obj: &mut SimobjT, env: &Environment) -> Array3d {
    let (r_0, r_tb, ob_0) = match sim_obj.state.soi {
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
            &sim_obj.state.coords_abs,
            &r_0,
            &r_tb.unwrap(),
            env.earth.attr.mass,
        )
    } else {
        newton_gravitational_field(&sim_obj.state.coords_abs, &r_0, env.earth.attr.mass)
    };

    let perturb = ecliptic_to_equatorial(&accel_ecliptic, ob_0);

    if let Some(perturb_store) = &mut sim_obj.perturb_store {
        perturb_store.earth_gravity = Some(PerturbationDefinition {
            perturb_name: "earth_gravity".to_string(),
            x_accel: perturb.x,
            y_accel: perturb.y,
            z_accel: perturb.z,
        });
    }

    perturb
}

pub fn calculate_earth_perturbations(sim_obj: &mut SimobjT, env: &Environment) -> Array3d {
    calculate_earth_gravity_perturbation(sim_obj, env)
        + calculate_earth_atmospheric_drag_perturbation(sim_obj, env)
}
