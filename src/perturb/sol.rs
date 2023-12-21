use crate::{
    bodies::{common::ecliptic_to_equatorial, sim_object::SimobjT, solar_model::Solarobj},
    environment::Environment,
    output,
    types::Array3d,
};

use super::common::{newton_gravitational_field, newton_gravitational_field_third_body};

fn calculate_solar_gravity_perturbation(
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    let (r_0, r_tb, ob_0) = match sim_obj.soi {
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
            &sim_obj.coords_abs,
            &r_0,
            &r_tb.unwrap(),
            env.sun.attr.mass,
        )
    } else {
        newton_gravitational_field(&sim_obj.coords_abs, &r_0, env.sun.attr.mass)
    };

    let accel = ecliptic_to_equatorial(&accel_ecliptic, ob_0);

    if let Some(out) = perturbations_out {
        out.push(output::PerturbationOut {
            id: sim_obj.id,
            sim_time: env.get_sim_time(),
            petrub_type: "solar_obj_sun".to_string(),
            acceleration_x_mpss: accel.x,
            acceleration_y_mpss: accel.y,
            acceleration_z_mpss: accel.z,
        })
    }

    accel
}

/// Calculate all solar perturbations on sim_obj
///
/// ### Arguments
/// * 'sim_obj' - Simulation object
/// * 'env' - Simulation environment
/// * 'perturbations_out' - Array of individual perturbations wrapped in a option.
///    This array should be appended to if option value is some
///
/// ### Returns
/// Combined perturbations
///
pub fn calculate_solar_perturbations(
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    calculate_solar_gravity_perturbation(sim_obj, env, perturbations_out)
}
