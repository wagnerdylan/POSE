use std::mem::swap;

use crate::bodies::sim_object::SimobjT;
use crate::collision::{collision_model, find_body_intersections, find_collision_set};
use crate::environment::Environment;
use crate::perturb;
use crate::{
    input,
    output::{self},
    types::Array3d,
};
use rayon::prelude::*;

/// Main function for calculating simulation object propagations.
///
/// ### Parameters
/// * 'sim_obj' - Simulation object in-which to apply propagations.
/// * 'env' - Simulation environment in which to perform propagations under.
/// * `step_time_s` - Simulation step-time value.
///
pub fn apply_perturbations(sim_obj: &mut SimobjT, env: &Environment, step_time_s: f64) {
    let mut net_acceleration = Array3d::default();
    // Calculate the perturbation forces for all planetary objects.
    net_acceleration = net_acceleration + perturb::sol::calculate_solar_perturbations(sim_obj, env);
    net_acceleration =
        net_acceleration + perturb::earth::calculate_earth_perturbations(sim_obj, env);
    net_acceleration = net_acceleration + perturb::moon::calculate_moon_perturbations(sim_obj, env);

    // Velocity and displacement calculations use Euler–Cromer integration as errors
    // in Euler–Cromer do not grow exponentially thus providing the most stable orbit.

    // Calculate the velocity change from net acceleration.
    let velocity_delta = net_acceleration * step_time_s;
    // Calculate new velocity for the given simulation object.
    let updated_sim_obj_velocity = velocity_delta + sim_obj.state.velocity;

    // Calculate the position change from the updated velocity.
    let position_displacement = updated_sim_obj_velocity * step_time_s;
    // Calculate the new position for the simulation object.
    let updated_sim_obj_coords = position_displacement + sim_obj.state.coords;

    // Update the new values within the simulation object.
    sim_obj.state.velocity = updated_sim_obj_velocity;
    sim_obj.state.coords = updated_sim_obj_coords;
}

/// Simulation halting condition check. Below are the conditions in which trigger simulation halting.
/// - Environment time equal to or exceeding configured stop time.
///
/// ### Parameters
/// * 'env' - Simulation environment used to perform propagation calculations.
/// * 'runtime_params' - Parameters collected on simulation init for running the sim.
///
/// ### Return
///     Boolean value indicating if the simulation should halt.
///
fn should_simulation_halt(env: &Environment, runtime_params: &input::RuntimeParameters) -> bool {
    env.start_time + chrono::Duration::milliseconds((env.sim_time_s * 1000.0) as i64)
        >= runtime_params.halt_date
}

/// Write out simulation results at the defined rate as configured.
///
/// ### Parameters
/// * 'sim_bodies' - Simulation objects to be propagated.
/// * 'env' - Simulation environment used to perform propagation calculations.
/// * 'output_controller' - Output object used for emitting simulation results.
/// * 'runtime_params' - Parameters collected on simulation init for running the sim.
/// * 'last_write' - Simulation time in which the last write occurred, -Infinity indicates this
///                  function should write results on call.
///
/// ### Return
///     On write of simulation data, return the current simulation time.
///     If data has not been written, return the 'last_write' value passed.
///
fn write_out_simulation_results(
    sim_bodies: &[SimobjT],
    env: &Environment,
    output_controller: &mut Box<dyn output::SimulationOutput>,
    runtime_params: &input::RuntimeParameters,
    last_write: f64,
) -> f64 {
    if env.sim_time_s < last_write + runtime_params.write_period - 0.001 {
        return last_write;
    }

    output::write_out_all_object_perturbations(env, sim_bodies, output_controller.as_mut());
    output::write_out_all_object_parameters(env, sim_bodies, output_controller.as_mut());
    output::write_out_all_solar_objects(env, output_controller.as_mut());

    // Return current sim time to allow for rated output.
    env.sim_time_s
}

fn propagate_simulation_objects(
    env: &Environment,
    runtime_params: &input::RuntimeParameters,
    sim_bodies: &mut [SimobjT],
) {
    // Calculate and apply perturbations for every object.
    sim_bodies.par_iter_mut().for_each(|sim_obj| {
        // Save previous solar ecliptic coordinates for intersection calculations.
        sim_obj.state.coord_helio_previous = sim_obj.state.coord_helio;
        // Update solar ecliptic coordinates for use in perturbation calculations.
        sim_obj.state.coord_helio = env.calculate_helio_coords(sim_obj);
        // Update fixed accelerating coordinates if applicable.
        sim_obj.state.coords_fixed = env.calculate_fixed_coords(sim_obj);

        env.check_switch_soi(sim_obj);
        apply_perturbations(sim_obj, env, runtime_params.sim_time_step as f64);
    });
}

fn run_collision_check(
    current_env: &Environment,
    previous_env: &mut Environment,
    runtime_params: &input::RuntimeParameters,
    sim_bodies: &mut Vec<SimobjT>,
) {
    // Step 01: find the set of overlapping bounding boxes which defines the
    // the collision set.

    // Return early if no collision sets exist.
    if !find_collision_set(sim_bodies) {
        return;
    }
    // Step 02: swap object parameters to that of the start of the collision intersection period.
    sim_bodies
        .iter_mut()
        .for_each(|a| swap(&mut a.state, &mut a.saved_state));

    while previous_env.step_count < current_env.step_count {
        let intersect_gen = sim_bodies
            .par_chunk_by_mut(|a, b| a.overlap_marker == b.overlap_marker)
            .flat_map(|slice| -> Vec<SimobjT> {
                // Skip over non-overlap group.
                if slice.first().unwrap().overlap_marker.is_none() {
                    return Vec::new();
                }
                // Step 03: propagate simulation objects by a single time step.
                propagate_simulation_objects(previous_env, runtime_params, slice);

                // Step 04: check intersections from previous simulation time step to most
                // recent time step propagation. If intersection results in collision, add
                // collision information into the list. sim_bodies may also be modified at this point
                // to reflect the result of the collision. IE. adding space debris into the simulation.
                let intersections = find_body_intersections(slice);
                intersections
                    .iter()
                    .flat_map(|intersect_idx| -> Vec<SimobjT> {
                        collision_model(
                            slice.get(intersect_idx.0).unwrap(),
                            slice.get(intersect_idx.1).unwrap(),
                        )
                    })
                    .collect()
            });
        // Append new generated simulation bodies into the main vector and re-sort to
        // chunk up simulation bodies by overlap_marker.
        let mut collision_gen: Vec<SimobjT> = intersect_gen.collect();
        if !collision_gen.is_empty() {
            sim_bodies.append(&mut collision_gen);
            sim_bodies.par_sort_unstable_by(|a, b| a.overlap_marker.cmp(&b.overlap_marker));
        }

        // End of simulation step calculations, prepare for next simulation step.
        previous_env.advance_simulation_environment(runtime_params);
    }

    sim_bodies
        .iter_mut()
        .for_each(|a| swap(&mut a.state, &mut a.saved_state));
}

fn should_run_collision_check(
    current_env: &Environment,
    previous_env: &Environment,
    runtime_params: &input::RuntimeParameters,
) -> bool {
    (current_env.sim_time_s - previous_env.sim_time_s)
        >= runtime_params.collision_check_period as f64
}

/// Main simulation loop
/// Exit from this loop is conditioned from simulation configuration.
///
/// ### Arguments
/// * 'sim_bodies' - Simulation objects to be propagated.
/// * 'env' - Initial simulation environment used to perform propagation calculations.
/// * 'output_controller' - Output object used for emitting simulation results.
/// * 'runtime_params' - Parameters collected on simulation init for running the sim.
///
/// ### Return
///     Return of this function is condition from simulation configuration.
///     This function does not return data.
///
pub fn simulate(
    mut sim_bodies: Vec<SimobjT>,
    mut env: Environment,
    mut output_controller: Box<dyn output::SimulationOutput>,
    runtime_params: input::RuntimeParameters,
) {
    let mut last_write = -f64::INFINITY;
    let mut previous_env = env.clone();

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

        // Main simulation object propagation.
        propagate_simulation_objects(&env, &runtime_params, &mut sim_bodies);

        // Check for simulation object collisions. Normal simulation state is restored after this
        // logic is run.
        if should_run_collision_check(&env, &previous_env, &runtime_params) {
            run_collision_check(&env, &mut previous_env, &runtime_params, &mut sim_bodies);
            previous_env = env.clone();
            sim_bodies
                .iter_mut()
                .for_each(|a| a.saved_state = a.state.clone());
        }

        // End of simulation step calculations, prepare for next simulation step.
        env.advance_simulation_environment(&runtime_params);

        // Simulation halting condition. Return from function "simulate" will end simulation execution.
        if should_simulation_halt(&env, &runtime_params) {
            return;
        }
    }
}
