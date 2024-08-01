use crate::{
    bodies::{
        common::{ecliptic_to_equatorial, spice_ephemeris_time},
        sim_object::{PerturbationDefinition, SimobjT},
        solar_model::Solarobj,
    },
    environment::Environment,
    types::{l2_norm, Array3d},
};

use chrono::{DateTime, Utc};
use satkit::{self, earthgravity::GravityModel};

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

/// High fidelity Earth gravity model which takes into account the J2
/// Earth gravity perturbation.
///
/// ### Arguments
/// * 'sim_obj_soi_coords' - Coordinates for a given simulation object expressed in the SOI.
/// * 'current_time' - Current simulation time object.
///
/// ### Return
///     Earth gravity acceleration for the given coordinates and time.
///
fn j2_gravity(sim_obj_soi_coords: &Array3d, current_time: &DateTime<Utc>) -> Array3d {
    let et = spice_ephemeris_time(current_time);
    // Rotate simulation object Earth coords in frame J2000 into ITRF93 for computing
    // zonal spherical harmonics Earth gravity.
    let j2000_irtf_rot_mtx = spice::pxform("J2000", "ITRF93", et);
    let pos_itrf = spice::mxv(j2000_irtf_rot_mtx, sim_obj_soi_coords.to_array());
    let j2_accel_itrf = satkit::earthgravity::accel(
        &satkit::jplephem::Vec3::new(
            *pos_itrf.get(0).unwrap(),
            *pos_itrf.get(1).unwrap(),
            *pos_itrf.get(2).unwrap(),
        ),
        8,
        GravityModel::JGM3,
    );
    // Accelerations returned by the spherical harmonics function are in the ITRF frame and
    // must be converted back into the J2000 frame.
    let irtf_j2000_rot_mtx = spice::pxform("ITRF93", "J2000", et);
    let j2_accel_itrf_array = [
        *j2_accel_itrf.get(0).unwrap(),
        *j2_accel_itrf.get(1).unwrap(),
        *j2_accel_itrf.get(2).unwrap(),
    ];
    let accel_j2000 = spice::mxv(irtf_j2000_rot_mtx, j2_accel_itrf_array);

    Array3d {
        x: *accel_j2000.get(0).unwrap(),
        y: *accel_j2000.get(1).unwrap(),
        z: *accel_j2000.get(2).unwrap(),
    }
}

fn calculate_earth_gravity_perturbation(sim_obj: &mut SimobjT, env: &Environment) -> Array3d {
    let (r_0, r_tb, ob_0) = match sim_obj.state.soi {
        Solarobj::Sun => (
            env.sun.model.state.coords,
            Some(env.earth.model.state.coords),
            env.sun.attr.obliquity,
        ),
        Solarobj::Earth => (env.earth.model.state.coords, None, env.earth.attr.obliquity),
        Solarobj::Moon => (
            env.moon.model.state.coords,
            Some(env.earth.model.state.coords),
            env.moon.attr.obliquity,
        ),
    };

    let perturb = if r_tb.is_some() {
        let accel = newton_gravitational_field_third_body(
            &sim_obj.state.coord_helio,
            &r_0,
            &r_tb.unwrap(),
            env.earth.attr.mass,
        );
        ecliptic_to_equatorial(&accel, ob_0)
    } else {
        // Use a lower fidelity but more preformat gravity model for debris.
        match sim_obj.sim_object_type {
            crate::bodies::sim_object::SimObjectType::Spacecraft => {
                j2_gravity(&sim_obj.state.coords, &env.current_time)
            }
            crate::bodies::sim_object::SimObjectType::Debris => {
                let accel = newton_gravitational_field(
                    &sim_obj.state.coord_helio,
                    &r_0,
                    env.earth.attr.mass,
                );
                ecliptic_to_equatorial(&accel, ob_0)
            }
        }
    };

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
