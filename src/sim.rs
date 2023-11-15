use crate::bodies::sim_object::SimobjT;
use crate::environment::Environment;
use crate::perturb;
use crate::{
    input,
    output::{self, PerturbationOut},
    types::Array3d,
};

const MAX_NUM_OF_PERTURBATIONS: usize = 4;

pub fn apply_perturbations(
    sim_obj: &mut SimobjT,
    env: &Environment,
    step_time_s: f64,
    perturbations_out: &mut Option<&mut Vec<PerturbationOut>>,
) {
    // Update absolute coordinates for use in perturbation calculations
    sim_obj.coords_abs = env.calculate_abs_coords(sim_obj);

    let mut net_acceleration = Array3d::default();
    // Calculate the perturbation forces for all planetary objects
    net_acceleration = net_acceleration
        + perturb::sol::calculate_solar_perturbations(sim_obj, env, perturbations_out);
    net_acceleration = net_acceleration
        + perturb::earth::calculate_earth_perturbations(sim_obj, env, perturbations_out);
    net_acceleration = net_acceleration
        + perturb::moon::calculate_moon_perturbations(sim_obj, env, perturbations_out);

    // Velocity and displacement calculations use Euler–Cromer integration as errors
    // in Euler–Cromer do not grow exponentially thus providing the most stable orbit.

    // Calculate the velocity change from net acceleration
    let velocity_delta = net_acceleration * step_time_s;
    // Calculate new velocity for the given simulation object
    let updated_sim_obj_velocity = velocity_delta + sim_obj.velocity;

    // Calculate the position change from the updated velocity
    let position_displacement = updated_sim_obj_velocity * step_time_s;
    // Calculate the new position for the simulation object
    let updated_sim_obj_coords = position_displacement + sim_obj.coords;

    // Update the new values within the simulation object
    sim_obj.velocity = updated_sim_obj_velocity;
    sim_obj.coords = updated_sim_obj_coords;
}

fn should_simulation_halt(env: &Environment, runtime_params: &input::RuntimeParameters) -> bool {
    env.start_time + chrono::Duration::milliseconds((env.get_sim_time() * 1000.0) as i64)
        >= runtime_params.halt_date
}

fn write_out_simulation_results(
    sim_bodies: &[SimobjT],
    env: &Environment,
    output_controller: &mut Box<dyn output::SimulationOutput>,
    runtime_params: &input::RuntimeParameters,
    last_write: f64,
) -> f64 {
    if env.get_sim_time() < last_write + runtime_params.write_period - 0.001 {
        return last_write;
    }

    output::write_out_all_object_parameters(env, sim_bodies, output_controller.as_mut());
    output::write_out_all_solar_objects(env, output_controller.as_mut());

    env.get_sim_time()
}

/// Main entry point into the cpu_sim module, gathers all needed data for orbit modeling
/// using Cowell's method.
pub fn simulate(
    mut sim_bodies: Vec<SimobjT>,
    mut env: Environment,
    mut output_controller: Box<dyn output::SimulationOutput>,
    runtime_params: input::RuntimeParameters,
) {
    // Allocate large buffer for holding perturbations
    let mut perturbation_vec: Option<Vec<PerturbationOut>> = match runtime_params.write_out_pertub {
        true => Some(Vec::with_capacity(
            sim_bodies.len() * MAX_NUM_OF_PERTURBATIONS,
        )),
        false => None,
    };
    let mut last_write = -f64::INFINITY;

    loop {
        // Write out simulation data at the start of the simulation loop as to capture
        // initial simulation state.
        last_write = write_out_simulation_results(
            &sim_bodies,
            &env,
            &mut output_controller,
            &runtime_params,
            last_write,
        );
        // Simulation halting condition. Return from function "simulate" will end simulation execution.
        if should_simulation_halt(&env, &runtime_params) {
            return;
        }

        // Calculate and apply perturbations for every object
        // TODO parallelize this
        for sim_obj in sim_bodies.iter_mut() {
            apply_perturbations(
                sim_obj,
                &env,
                runtime_params.sim_time_step as f64,
                &mut perturbation_vec.as_mut(),
            );

            // TODO only call this every so often
            env.check_switch_soi(sim_obj);
        }

        if let Some(vec) = &mut perturbation_vec {
            output::write_out_all_perturbations(vec, output_controller.as_mut());
            vec.clear();
        }

        env.advance_simulation_environment(&runtime_params);
    }
}
