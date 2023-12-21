use crate::{
    bodies::{common::ecliptic_to_equatorial, sim_object::SimobjT, solar_model::Solarobj},
    environment::Environment,
    output,
    types::Array3d,
};

use super::common::{newton_gravitational_field, newton_gravitational_field_third_body};

fn calculate_moon_gravity_perturbation(
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
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

    let accel = ecliptic_to_equatorial(&accel_ecliptic, ob_0);

    if let Some(out) = perturbations_out {
        out.push(output::PerturbationOut {
            id: sim_obj.id,
            sim_time: env.get_sim_time(),
            petrub_type: "solar_obj_moon".to_string(),
            acceleration_x_mpss: accel.x,
            acceleration_y_mpss: accel.y,
            acceleration_z_mpss: accel.z,
        })
    }

    accel
}

pub fn calculate_moon_perturbations(
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    calculate_moon_gravity_perturbation(sim_obj, env, perturbations_out)
}
