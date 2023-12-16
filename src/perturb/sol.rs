use crate::{
    bodies::{
        self,
        sim_object::SimobjT,
        solar_model::{KeplerModel, Solarobj},
    },
    environment::Environment,
    output,
    types::Array3d,
};
use perturb::common;

fn calculate_solar_gravity_perturbation(
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    // Calculate distance between the sun and the sim object using absolute coordinates
    let distance_vector_sim_obj = sim_obj.coords_abs - env.sun.model.state.coords.current_coords;

    // Calculate acceleration due to the sun at the location of the simulation object
    let mut gravity_accel = common::newton_gravitational_field(
        &distance_vector_sim_obj,
        env.sun.model.get_solar_object().get_mass_kg(),
    );

    // If the simulation object is not in the solar SOI, the gravity field must be compensated
    // This is because perturbation calculations are done relative to the SOI
    gravity_accel = match sim_obj.soi {
        Solarobj::Sun { attr: _ } => gravity_accel,
        Solarobj::Earth { attr: _ } => bodies::common::ecliptic_to_equatorial(
            &(gravity_accel
                - common::newton_gravitational_field(
                    &(env.earth.model.state.coords.current_coords
                        - env.sun.model.state.coords.current_coords),
                    env.sun.model.get_solar_object().get_mass_kg(),
                )),
            env.earth.model.get_solar_object().get_obliquity(),
        ),
        Solarobj::Moon { attr: _ } => bodies::common::ecliptic_to_equatorial(
            &(gravity_accel
                - common::newton_gravitational_field(
                    &(env.moon.model.state.coords.current_coords
                        - env.sun.model.state.coords.current_coords),
                    env.sun.model.get_solar_object().get_mass_kg(),
                )),
            env.moon.model.get_solar_object().get_obliquity(),
        ),
    };

    if let Some(out) = perturbations_out {
        out.push(output::PerturbationOut {
            id: sim_obj.id,
            sim_time: env.get_sim_time(),
            petrub_type: "solar_obj_sun".to_string(),
            acceleration_x_mpss: gravity_accel.x,
            acceleration_y_mpss: gravity_accel.y,
            acceleration_z_mpss: gravity_accel.z,
        })
    }

    gravity_accel
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
