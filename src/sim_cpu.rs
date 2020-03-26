use crate::bodies;

struct PerturbationDelta {
    acceleration_x_mpss: f32,
    acceleration_y_mpss: f32,
    acceleration_z_mpss: f32,
}

impl Default for PerturbationDelta {
    fn default() -> Self {
        Self {
            acceleration_x_mpss: 0.0,
            acceleration_y_mpss: 0.0,
            acceleration_z_mpss: 0.0,
        }
    }
}

pub enum Perturbation {
    SolarObject(bodies::Solarobj, PerturbationDelta),
}

/// Module used to apply perturbation calculations on individual bodies
mod cowell_perturb {
    use crate::bodies;
    use crate::sim_cpu::Perturbation;

    /// Apply all perturbations handled by POSE. This includes:
    /// * 'Solar Body Earth'
    /// * 'Solar Body Moon'
    /// * 'Solar Body Sun'
    /// TODO add more
    ///
    /// ### Parameters
    /// * 'sim_obj' - The object basis for calculation and apply
    /// * 'env' - The Simulation environment
    ///
    /// ### Return
    ///     A vector of all perturbations applied to the object basis.
    ///
    pub fn apply_perturbations(
        sim_obj: &mut bodies::SimobjT,
        env: &bodies::Environment,
    ) -> Vec<Perturbation> {
        // TODO collect a vector of all perturbations
        // TODO sum perturbations
        // TODO apply sum of perturbations to sim_obj
        // TODO emit vector of all perturbations

        unimplemented!();
    }

    /// Calculate perturbations due to solar system objects.
    ///
    /// ### Parameters
    /// * 'sim_obj' - The object basis for calculation
    /// * 'env' - The Simulation environment
    ///
    /// ### Return
    ///     A vector of all solar system perturbations calculated for this object.
    ///
    fn calc_planet_perturb(
        sim_obj: &bodies::SimobjT,
        env: &bodies::Environment,
    ) -> Vec<Perturbation> {
        // TODO
        unimplemented!();
    }

}

/// Main entry point into the cpu_sim module, gathers all needed data for orbit modeling
/// using Cowell's method.
pub fn simulate(mut sim_bodies: Vec<bodies::SimobjT>, mut env: bodies::Environment) {}
