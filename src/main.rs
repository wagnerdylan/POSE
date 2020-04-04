//!
//! POSE - Parallel Orbital Simulation Environment
//! TODO - Add more doc

#[macro_use]
mod macros;

extern crate chrono;
extern crate clap;
extern crate serde;
extern crate serde_json;
extern crate strum;
extern crate strum_macros;
extern crate csv;

mod bodies;
mod input;
mod output;
mod sim_cpu;

mod cli {

    ///Checks if value passed in to program argument is numeric. Returns a Result
    ///
    ///# Argument
    ///* 'strng' - The value passed by the user
    ///
    fn numeric_validator(strng: String) -> Result<(), String> {
        if strng.parse::<f32>().is_ok() {
            Ok(())
        } else {
            Err(String::from("Input is non-numeric"))
        }
    }

    /// Defines the argument structure for the pose simulation program
    /// Returns the result of user arguments passed over the cli
    pub fn check_cli() -> clap::ArgMatches<'static> {
        // Defines the input arguments from the cli
        let matches = clap::App::new("Parallel Orbital Simulation Environment (POSE)")
            .version("DEV0.1")
            .about("Simulation aimed to model the orbital environment around Earth for bodies at all magnitudes.")
            .args(&[
                clap::Arg::with_name("INPUT")
                    .help("json file containing information on bodies at initialization.")
                    .required(true)
                    .index(1),
                clap::Arg::with_name("out")
                    .help("Output specifier")
                    .short("o")
                    .long("out")
                    .value_name("DIR_NAME")
                    .takes_value(true),
                clap::Arg::with_name("sim_time_step")
                    .help("Simulation time step interval in seconds")
                    .short("s")
                    .long("step")
                    .value_name("STEP_INTERVAL")
                    .takes_value(true)
                    .validator(numeric_validator)
            ])
            .get_matches();

        return matches;
    }
}

fn main() {
    let matches = cli::check_cli();
    let sim_params = input::gather_program_arguments(matches);

    let (sim_bodies, start_time) =
        input::parse_input(sim_params.input_bodies_json.as_str());
    let env = bodies::Environment::new(start_time);

    let output_controller =
        Box::new(output::csv_output::CSVController::new(sim_params.output_dir.as_str()));

    sim_cpu::simulate(sim_bodies, env, output_controller);
}
