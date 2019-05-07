use crate::bodies;

/// Module used to apply perturbation calculations on and individual body
mod cowell_perturb {
    use crate::bodies;

    fn planet_perturb(sim_obj: bodies::SimobjT, solar_bodies: &Vec<bodies::PlanetBody>)
        -> bodies::SimobjT
    {
        // TODO
        sim_obj
    }

}

/// Main entry point into the cpu_sim module, gathers all needed data for orbit modeling
/// using Cowell's method.
pub fn simulate(
    sim_bodies: Vec<bodies::SimobjT>,
    env: bodies::Environment,
    day: f32,
    step: f32,
    ){


}