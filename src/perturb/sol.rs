use crate::{
    bodies::{
        common::ecliptic_to_equatorial,
        sim_object::{PerturbationDefinition, SimobjT},
        solar_model::Solarobj,
    },
    environment::Environment,
    types::Array3d,
};

use super::common::{newton_gravitational_field, newton_gravitational_field_third_body};

fn calculate_solar_gravity_perturbation(sim_obj: &mut SimobjT, env: &Environment) -> Array3d {
    let (r_0, r_tb, ob_0) = match sim_obj.state.soi {
        Solarobj::Sun => (
            env.sun.model.state.coords.current_coords,
            None,
            env.sun.attr.obliquity,
        ),
        Solarobj::Earth => (
            env.earth.model.state.coords.current_coords,
            Some(env.sun.model.state.coords.current_coords),
            env.earth.attr.obliquity,
        ),
        Solarobj::Moon => (
            env.moon.model.state.coords.current_coords,
            Some(env.sun.model.state.coords.current_coords),
            env.moon.attr.obliquity,
        ),
    };

    let accel_ecliptic = if r_tb.is_some() {
        newton_gravitational_field_third_body(
            &sim_obj.state.coords_abs,
            &r_0,
            &r_tb.unwrap(),
            env.sun.attr.mass,
        )
    } else {
        newton_gravitational_field(&sim_obj.state.coords_abs, &r_0, env.sun.attr.mass)
    };

    let perturb = ecliptic_to_equatorial(&accel_ecliptic, ob_0);

    if let Some(perturb_store) = &mut sim_obj.perturb_store {
        perturb_store.sun_gravity = Some(PerturbationDefinition {
            perturb_name: "sun_gravity".to_string(),
            x_accel: perturb.x,
            y_accel: perturb.y,
            z_accel: perturb.z,
        });
    }

    perturb
}

/// Calculate all solar perturbations on sim_obj
///
/// ### Arguments
/// * 'sim_obj' - Simulation object
/// * 'env' - Simulation environment
/// ### Returns
/// Combined perturbations
///
pub fn calculate_solar_perturbations(sim_obj: &mut SimobjT, env: &Environment) -> Array3d {
    calculate_solar_gravity_perturbation(sim_obj, env)
}
