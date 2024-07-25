use std::mem::swap;

use crate::bodies::sim_object::SimobjT;
use crate::collision::{
    collision_model, find_body_intersections, find_collision_set, CollisionResult,
};
use crate::environment::Environment;
use crate::perturb;
use crate::types::l2_norm;
use crate::{
    input,
    output::{self},
    types::Array3d,
};
use rayon::prelude::*;
use std::slice;

/// Mark a given simulation object for deletion if the object has dipped
/// below the surface of the SOI body. Objects marked by this function may be
/// culled later on.
///
/// ### Arguments
/// * 'sim_obj' - Simulation object to be checked for deletion.
/// * 'env' - Simulation environment which matches the point in time in which sim_obj exists.
///
pub fn check_soi_intersection(sim_obj: &mut SimobjT, env: &Environment) {
    let soi_radius = {
        match sim_obj.state.soi {
            crate::bodies::solar_model::Solarobj::Sun => env.sun.attr.eqradius,
            crate::bodies::solar_model::Solarobj::Earth => env.earth.attr.eqradius,
            crate::bodies::solar_model::Solarobj::Moon => env.moon.attr.eqradius,
        }
    };

    let distance_to_soi_body = l2_norm(&(sim_obj.state.coords));
    if distance_to_soi_body < soi_radius {
        sim_obj.marked_for_deletion_on_step = Some(env.step_count);
    }
}

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

fn update_derived_coords(sim_obj: &mut SimobjT, env: &Environment) {
    sim_obj.state.coord_helio = env.calculate_helio_coords(sim_obj);
    sim_obj.state.coords_fixed = env.calculate_fixed_coords(sim_obj);
}

/// Main propagation function for simulation objects.
/// This function modifies state of each simulation object in-place.
///
/// ### Parameters
/// * 'env' - Environment to use in object propagation.
/// * 'runtime_params' - Parameters collected on simulation init for running the sim.
/// * 'sim_bodies' - Simulation objects to be propagated.
///
fn propagate_simulation_objects(
    env: &Environment,
    runtime_params: &input::RuntimeParameters,
    sim_bodies: &mut [SimobjT],
) {
    // Calculate and apply perturbations for every object.
    sim_bodies.par_iter_mut().for_each(|sim_obj| {
        // Save previous solar ecliptic coordinates for intersection calculations.
        sim_obj.state.previous_coords = sim_obj.state.coords;

        // Update derived coordinates to decouple env update call order.
        update_derived_coords(sim_obj, env);

        env.check_switch_soi(sim_obj);
        check_soi_intersection(sim_obj, env);
        apply_perturbations(sim_obj, env, runtime_params.sim_time_step as f64);

        // Update derived coordinates after object propagations have been applied.
        // After this update all object coordinates will be in-sync.
        update_derived_coords(sim_obj, env);
    });
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

/// Main collision check algorithm.
///
/// (1) Find overlapping bounding boxes which define a group of objects
///     to be checked for intersections. Groups are layed out within the
///     sim_bodies vector sequentially. For instance sim_bodies may look
///     like [None, None, 1, 1, 2, 2, 2...]. None markers indicate that
///     a given simulation object does not belong to a intersection group.
///     The 'None' group is ignored in the following steps.
/// (2) Prepare for object re-propagation by swapping the saved_state with
///     the current state for objects which belong to a collision group.
///     This effectively 'rewinds' the simulation so intersections may be checked.
/// (3) Propagate simulation objects using the previous_env checking for
///     intersections between each object within a intersection group.
///     If an intersection occurs, the collision model is called which
///     may generate new simulation bodies and delete current objects
///     passed into the model.
///
/// ### Pre-Conditions
///     'current_env' is at the bound of the most recent collision check period.
///     'previous_env' is at the bound of the last collision check period
///     (this may also be the starting point of the simulation).
///     'sim_bodies' has perturbation calculations applied up to the point of the
///     most recent collision check period (at the same point as current_env).
///
/// ### Post-Conditions
///     'previous_env' may or may not have been updated to the point of 'current_env'
///     as such, the state of 'previous_env' should not be relied upon return.
///     Simulation objects within 'sim_bodies' which are marked for collision checking
///     will be re-propagated until the point of 'saved_state'.
///
/// ### Parameters
/// * 'current_env' - Most recent simulation environment used as a stopping point.
/// * 'previous_env' - Simulation environment from the last collision check period.
/// * 'runtime_params' - Parameters collected on simulation init for running the sim.
/// * 'sim_bodies' - Simulation objects to be propagated.
/// * 'output_controller' - Output object used for emitting simulation results.
///
fn run_collision_check(
    current_env: &Environment,
    previous_env: &mut Environment,
    runtime_params: &input::RuntimeParameters,
    sim_bodies: &mut Vec<SimobjT>,
    output_controller: &mut Box<dyn output::SimulationOutput>,
) {
    // Return early if no collision sets exist.
    if !find_collision_set(sim_bodies, runtime_params.check_only_satellite_collisions) {
        return;
    }

    let max_id = sim_bodies.iter().max_by(|a, b| a.id.cmp(&b.id)).unwrap().id;
    let mut new_max_id = max_id;

    // Swap in the saved state which is from the beginning of the main simulation collision check period.
    // This is done only to the simulation objects which will be checked for intersections.
    sim_bodies
        .chunk_by_mut(|a, b| a.overlap_marker == b.overlap_marker)
        .for_each(|slice| {
            if !slice.first().unwrap().overlap_marker.is_none() {
                slice
                    .iter_mut()
                    .for_each(|a| swap(&mut a.state, &mut a.saved_state))
            }
        });

    while previous_env.step_count < current_env.step_count {
        let intersect_gen = sim_bodies
            .par_chunk_by_mut(|a, b| a.overlap_marker == b.overlap_marker)
            .flat_map(|slice| -> Vec<CollisionResult> {
                // Skip over non-overlap group.
                if slice.first().unwrap().overlap_marker.is_none() {
                    return Vec::new();
                }

                // Single time-step object propagation.
                propagate_simulation_objects(previous_env, runtime_params, slice);

                // Check intersections from previous simulation time step to most
                // recent time step propagation. If intersection results in collision, add
                // collision information into the list. sim_bodies may also be modified at this point
                // to reflect the result of the collision. IE. adding space debris into the simulation.
                let intersections = find_body_intersections(slice, previous_env.step_count);
                intersections
                    .iter()
                    .map(|intersect_result| -> CollisionResult {
                        collision_model(previous_env, slice, intersect_result)
                    })
                    .collect()
            });

        let collision_gen: Vec<CollisionResult> = intersect_gen.collect();
        if !collision_gen.is_empty() {
            // Include newly generated simulation objects into the main sim_bodies array.
            let mut new_bodies: Vec<SimobjT> = collision_gen
                .iter()
                .flat_map(|a| a.new_sim_bodies.to_owned())
                .collect();
            // Id values must not be duplicated and ascend without sequence gaps.
            new_bodies.iter_mut().for_each(|a| {
                new_max_id += 1;
                a.id = new_max_id
            });
            sim_bodies.append(&mut new_bodies);
            // Sort by overlap marker here so generated bodies are placed within the objects which
            // collided to generated them.
            sim_bodies.par_sort_unstable_by(|a, b| a.overlap_marker.cmp(&b.overlap_marker));

            collision_gen.iter().for_each(|x| {
                output_controller.write_out_collision_info(x.to_output_form(&new_bodies))
            })
        }

        // Write out simulation results for new objects to populate simulation
        // output so no gaps in information occur. Write output at the simulation
        // step rate.
        if new_max_id != max_id {
            for sim_body in sim_bodies.iter() {
                if sim_body.id > max_id {
                    write_out_simulation_results(
                        slice::from_ref(sim_body),
                        previous_env,
                        output_controller,
                        runtime_params,
                        -f64::INFINITY,
                    );
                }
            }
        }

        // End of simulation step calculations, prepare for next simulation step.
        previous_env.advance_simulation_environment(runtime_params);
    }

    // At the end of the function, saved_state and state within each object checked for
    // collision should be the same.
}

fn should_run_sim_reverse_logic(
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
        // Simulation halting condition. Return from function "simulate" will end simulation execution.
        if should_simulation_halt(&env, &runtime_params) {
            return;
        }

        // Write out simulation data at the start of the simulation loop as to capture
        // initial simulation state.
        last_write = write_out_simulation_results(
            &sim_bodies,
            &env,
            &mut output_controller,
            &runtime_params,
            last_write,
        );

        // Check for simulation object collisions. Normal simulation state is restored after this
        // logic is run.
        if should_run_sim_reverse_logic(&env, &previous_env, &runtime_params) {
            run_collision_check(
                &env,
                &mut previous_env,
                &runtime_params,
                &mut sim_bodies,
                &mut output_controller,
            );
            // Remove all bodies marked for deletion at this point. The step count contained within
            // the marked_for_deletion_on_step member is irreverent at this point as it is always in the past.
            sim_bodies.retain(|a| a.marked_for_deletion_on_step.is_none());
            previous_env = env.clone();
            sim_bodies
                .iter_mut()
                .for_each(|a| a.saved_state = a.state.clone());
        }

        // Main simulation object propagation.
        propagate_simulation_objects(&env, &runtime_params, &mut sim_bodies);

        // End of simulation step calculations, prepare for next simulation step.
        env.advance_simulation_environment(&runtime_params);
    }
}
