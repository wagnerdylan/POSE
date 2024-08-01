//!
//! POSE - Parallel Orbital Simulation Environment
//!

use bodies::{sim_object::SimobjT, solar_model::furnish_spice};
use clap::Parser;

#[macro_use]
mod macros;
#[macro_use]
extern crate impl_ops;

extern crate chrono;
extern crate clap;
extern crate csv;
extern crate nalgebra;
extern crate nrlmsise00c;
extern crate rayon;
extern crate satkit;
extern crate serde;
extern crate serde_json;
extern crate spice;
extern crate strum;
extern crate strum_macros;

mod bodies;
mod collision;
mod environment;
mod input;
mod output;
mod perturb;
mod sim;
mod types;

/// Initialize derived information from source information for each object within the simulation.
///
/// ### Arguments
/// * 'env' - Simulation environment
/// * 'sim_objs' - Array slice of simulation objects
///
fn init_simulation_objects(env: &environment::Environment, sim_objs: &mut Vec<SimobjT>) {
    for sim_obj in sim_objs {
        sim_obj.state.coord_helio = env.calculate_helio_coords(sim_obj);
        sim_obj.state.coords_fixed = env.calculate_fixed_coords(sim_obj);
        sim_obj.saved_state = sim_obj.state.clone();
    }
}

fn main() {
    let sim_params = input::SimulationParameters::parse();
    // furnish_spice must be called before initializing solar objects.
    furnish_spice();

    let (mut sim_obj_holder, runtime_params, env_init) =
        input::collect_simulation_inputs(&sim_params);
    let env = environment::Environment::new(&runtime_params, env_init);

    init_simulation_objects(&env, &mut sim_obj_holder.sim_objs);

    let output_controller = Box::new(output::csv_output::CSVController::new(
        sim_params.output_dir.as_str(),
    ));

    sim::simulate(sim_obj_holder, env, output_controller, runtime_params);
}
