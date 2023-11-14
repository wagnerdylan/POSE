use crate::types::{self, Array3d};

// Gravitational constant 6.674×10−11
const G: f64 = 6.674e-11;
pub fn newton_gravitational_field(distance_vector: &Array3d, planet_mass_kg: f64) -> Array3d {
    let l2_dist = types::l2_norm(distance_vector);
    // Calculate unit vector for perturbation
    let unit_vector = types::normalize(distance_vector, Some(l2_dist));
    // Calculate acceleration field using Newton's law of universal gravitation
    unit_vector * (-G * (planet_mass_kg / l2_dist.powi(2)))
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
            1e15f64,
        );

        assert!(-12845.0 < grav_result.x && grav_result.x < -12843.0);
        assert!(-12845.0 < grav_result.y && grav_result.y < -12843.0);
        assert!(-12845.0 < grav_result.z && grav_result.z < -12843.0);
    }
}
