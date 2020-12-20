use crate::bodies;

use super::sim_cpu;

fn calculate_solar_gravity_perturbation(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
) -> sim_cpu::Perturbation {
    todo!()
}

/// Calculate all solar perturbations on sim_obj
///
/// ### Arguments
/// * 'sim_obj' - Simulation object
/// * 'env' - Simulation environment
/// * 'perturbations_out' - Array of individual perturbations wrapped in a option.
///    This array should be appended to if option value is some
/// * 'step_time_s' - Simulation step time
///
/// ### Returns
/// Combined perturbations
///
pub fn calculate_solar_perturbations(
    sim_obj: &bodies::SimobjT,
    env: &bodies::Environment,
    perturbations_out: Option<&mut [sim_cpu::Perturbation]>,
    step_time_s: f64,
) -> sim_cpu::PerturbationDelta {
    todo!()
}
