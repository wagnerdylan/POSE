//!
//! POSE - Parallel Orbital Simulation Environment
//! TODO - Add more doc

use bodies::{Environment, SimobjT};
use clap::Parser;

#[macro_use]
mod macros;
#[macro_use]
extern crate impl_ops;

extern crate chrono;
extern crate clap;
extern crate csv;
extern crate nrlmsise00c;
extern crate serde;
extern crate serde_json;
extern crate strum;
extern crate strum_macros;

mod bodies;
mod cpu;
mod input;
mod output;
mod types;

/// Initialize derived information from source information for each object within the simulation.
///
/// ### Arguments
/// * 'env' - Simulation environment
/// * 'sim_objs' - Array slice of simulation objects
///
fn init_simulation_objects(env: &Environment, sim_objs: &mut Vec<SimobjT>) {
    for sim_obj in sim_objs {
        sim_obj.coords_abs = env.calculate_abs_coords(sim_obj);
    }
}

fn main() {
    let sim_params = input::SimulationParameters::parse();

    let (mut sim_bodies, start_time) = input::parse_input(sim_params.configuration.as_str());
    let env = bodies::Environment::new(start_time, &sim_params);

    init_simulation_objects(&env, &mut sim_bodies);

    // TODO: dont box this
    let output_controller = Box::new(output::csv_output::CSVController::new(
        sim_params.output_dir.as_str(),
    ));

    cpu::sim_cpu::simulate(sim_bodies, env, output_controller, sim_params);
}
