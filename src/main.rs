//!
//! POSE - Parallel Orbital Simulation Environment
//! TODO - Add more doc

#[macro_use]
mod macros;

extern crate clap;

mod innout;
mod bodies;
mod sim_cpu;

mod cli{

    ///Checks if value passed in to program argument is numeric. Returns a Result
    ///
    ///# Argument
    ///* 'strng' - The value passed by the user
    fn numeric_validator(strng: String) -> Result<(), String>{
        if strng.parse::<f32>().is_ok(){
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
                    .help("json file containing information on bodies at initilzation.")
                    .required(true)
                    .index(1),
                clap::Arg::with_name("out")
                    .help("Main name of output files.")
                    .short("o")
                    .long("out")
                    .value_name("FILE_NAME")
                    .takes_value(true),
                clap::Arg::with_name("step")
                    .help("Simulation time step interval in seconds")
                    .short("s")
                    .long("step")
                    .value_name("STEP_INTERVAL")
                    .takes_value(true)
                    .default_value("0.1") // Simulation step time
                    .validator(numeric_validator)
            ])
            .get_matches();

        return matches;
    }
}

fn main() {
    let matches = cli::check_cli();
    let inpt_file = matches.value_of("INPUT").unwrap(); // Will always have INPUT
    let step = matches.value_of("step").unwrap().parse::<f32>().unwrap(); // Always have default

    let (sim_bodies, day) = innout::parse_inpt(inpt_file);
    let env = bodies::solar_system_objs(day);

    sim_cpu::simulate(sim_bodies, env, day, step);
}
