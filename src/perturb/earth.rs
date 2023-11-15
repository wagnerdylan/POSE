use crate::{
    bodies::{
        sim_object::SimobjT,
        solar_object::{KeplerModel, Solarobj},
    },
    environment::Environment,
    output,
    types::{l2_norm, Array3d},
};
use chrono::{Datelike, Duration, Timelike};
use perturb::common;

fn nrlmsise00_model(env: &Environment, sim_obj_alt_meters: f64) -> nrlmsise00c::NRLMSISEOutput {
    let current_datetime =
        env.start_time + Duration::milliseconds((env.get_sim_time() * 1000.0) as i64);
    let seconds_from_midnight = current_datetime.num_seconds_from_midnight() as f64;
    let mut input = nrlmsise00c::NRLMSISEInput {
        year: current_datetime.year(),
        doy: current_datetime.ordinal() as i32,
        sec: seconds_from_midnight,
        alt: sim_obj_alt_meters / 1000.0,
        // Use location of prime meridian as static point.
        // This may be enhanced by converting sim obj from ECI to ECEF.
        g_lat: 0.0,
        g_long: 0.0,
        lst: seconds_from_midnight / 3600.0 + 0.0 / 15.0, // (lst=sec/3600 + g_long/15)
        f107A: 69.0,
        f107: 69.0,
        ap: 7.0,
        ap_a: [0f64; 7usize],
    };
    let flags = nrlmsise00c::NRLMSISEFlags {
        switches: [1; 24usize],
        sw: [0f64; 24usize],
        swc: [0f64; 24usize],
    };

    nrlmsise00c::gtd7_safe(&mut input, &flags)
}

fn calculate_earth_atmospheric_drag_perturbation(
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    let distance_sim_obj = l2_norm(&(sim_obj.coords_abs - env.earth.coords.current_coords));
    let sim_obj_alt = distance_sim_obj - env.earth.get_solar_object().get_radius_meters();

    // return an empty result if sim object earth relative altitude is not within a range which
    // atmospheric drag will be relevant to simulation.
    if !(0.0..=20_000.0 * 1000.0).contains(&sim_obj_alt) {
        return Array3d::default();
    }

    let atmos_model_result = nrlmsise00_model(env, sim_obj_alt);
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
            petrub_type: "drag_earth".to_string(),
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
    let distance_vector_sim_obj = sim_obj.coords_abs - env.earth.coords.current_coords;

    // Calculate acceleration due to the earth at the location of the simulation object
    let mut gravity_accel = common::newton_gravitational_field(
        &distance_vector_sim_obj,
        env.earth.get_solar_object().get_mass_kg(),
    );

    // If the simulation object is not in the solar SOI, the gravity field must be compensated
    // This is because perturbation calculations are done relative to the SOI
    let soi_compensation = match sim_obj.soi {
        Solarobj::Sun { attr: _ } => None,
        Solarobj::Earth { attr: _ } => None,
        Solarobj::Moon { attr: _ } => {
            Some(env.moon.coords.current_coords - env.earth.coords.current_coords)
        }
    };

    // Calculate the difference between the field on simulation object and solar object
    if let Some(soi_info) = soi_compensation {
        gravity_accel = gravity_accel
            - common::newton_gravitational_field(
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
    sim_obj: &SimobjT,
    env: &Environment,
    perturbations_out: &mut Option<&mut Vec<output::PerturbationOut>>,
) -> Array3d {
    calculate_earth_gravity_perturbation(sim_obj, env, perturbations_out)
        + calculate_earth_atmospheric_drag_perturbation(sim_obj, env, perturbations_out)
}
