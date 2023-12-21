use crate::types::{self, Array3d};

// Gravitational constant 6.674×10−11
const G: f64 = 6.674e-11;

/// Point mass gravity field without third body.
///
/// ### Arguments
/// * 'sim_obj_coord' - Coordinates of simulation object.
/// * 'primary_body_coord' - Coordinates of the primary body.
/// * 'planet_mass_kg' - Mass of the primary body.
///
/// ### Returns
///     Acceleration field of the primary body.
///
pub fn newton_gravitational_field(
    sim_obj_coord: &Array3d,
    primary_body_coord: &Array3d,
    planet_mass_kg: f64,
) -> Array3d {
    let r = sim_obj_coord - primary_body_coord;
    let u = G * planet_mass_kg;
    // Calculate acceleration field using Newton's law of universal gravitation
    -u * (r / types::l2_norm(&r).powi(3))
}

/// Point mass gravity field with third body compensation.
/// Ref: https://www.agi.com/resources/whitepapers/correct-modeling-of-the-indirect-term-for-third-bo
///
/// ### Arguments
/// * 'sim_obj_coord' - Coordinates of simulation object.
/// * 'primary_body_coord' - Coordinates of the primary body.
/// * 'third_body_coord' - Coordinates of third body.
/// * 'planet_mass_kg' - Mass of the primary body.
///
/// ### Returns
///     Compensated gravity field for third body.
///
pub fn newton_gravitational_field_third_body(
    sim_obj_coord: &Array3d,
    primary_body_coord: &Array3d,
    third_body_coord: &Array3d,
    planet_mass_kg: f64,
) -> Array3d {
    let r = sim_obj_coord - primary_body_coord;
    let r_b = third_body_coord - primary_body_coord;
    let u = G * planet_mass_kg;

    let r1 = r_b - r;
    u * (r1 / types::l2_norm(&r1).powi(3) - r_b / types::l2_norm(&r_b).powi(3))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_newton_gravitational_field() {
        let grav_result = newton_gravitational_field(
            &Array3d {
                x: 1.0,
                z: 1.0,
                y: 1.0,
            },
            &Array3d {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            1e15f64,
        );

        assert!(-12845.0 < grav_result.x && grav_result.x < -12843.0);
        assert!(-12845.0 < grav_result.y && grav_result.y < -12843.0);
        assert!(-12845.0 < grav_result.z && grav_result.z < -12843.0);
    }
}
