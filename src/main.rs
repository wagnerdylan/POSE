//!
//! POSE - Parallel Orbital Simulation Environment
//!

use bodies::sim_object::SimobjT;
use clap::Parser;

#[macro_use]
mod macros;
#[macro_use]
extern crate impl_ops;

extern crate chrono;
extern crate clap;
extern crate csv;
extern crate nrlmsise00c;
extern crate rayon;
extern crate serde;
extern crate serde_json;
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
        sim_obj.state.coord_helio = env.calculate_se_coords(sim_obj);
        sim_obj.state.coords_fixed = env.calculate_fixed_coords(sim_obj);
        sim_obj.saved_state = sim_obj.state.clone();
    }
}

fn main() {
    let sim_params = input::SimulationParameters::parse();

    let (mut sim_bodies, runtime_params, env_init) = input::collect_simulation_inputs(&sim_params);
    let env = environment::Environment::new(&runtime_params, env_init);

    init_simulation_objects(&env, &mut sim_bodies);

    let output_controller = Box::new(output::csv_output::CSVController::new(
        sim_params.output_dir.as_str(),
    ));

    sim::simulate(sim_bodies, env, output_controller, runtime_params);
}
