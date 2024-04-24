use crate::{
    bodies::{common::ecliptic_to_equatorial, sim_object::SimobjT, solar_model::Solarobj},
    environment::Environment,
    types::Array3d,
};

use super::common::{newton_gravitational_field, newton_gravitational_field_third_body};

fn calculate_moon_gravity_perturbation(sim_obj: &SimobjT, env: &Environment) -> Array3d {
    let (r_0, r_tb, ob_0) = match sim_obj.soi {
        Solarobj::Sun => (
            env.sun.model.state.coords.current_coords,
            Some(env.moon.model.state.coords.current_coords),
            env.sun.attr.obliquity,
        ),
        Solarobj::Earth => (
            env.earth.model.state.coords.current_coords,
            Some(env.moon.model.state.coords.current_coords),
            env.earth.attr.obliquity,
        ),
        Solarobj::Moon => (
            env.moon.model.state.coords.current_coords,
            None,
            env.moon.attr.obliquity,
        ),
    };

    let accel_ecliptic = if r_tb.is_some() {
        newton_gravitational_field_third_body(
            &sim_obj.coords_abs,
            &r_0,
            &r_tb.unwrap(),
            env.moon.attr.mass,
        )
    } else {
        newton_gravitational_field(&sim_obj.coords_abs, &r_0, env.moon.attr.mass)
    };

    ecliptic_to_equatorial(&accel_ecliptic, ob_0)
}

pub fn calculate_moon_perturbations(sim_obj: &SimobjT, env: &Environment) -> Array3d {
    calculate_moon_gravity_perturbation(sim_obj, env)
}
