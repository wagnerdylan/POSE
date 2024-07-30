use chrono::{DateTime, Duration, Utc};

use crate::{
    bodies::{
        self,
        common::{days_since_j2000, spice_ephemeris_time},
        sim_object::SimobjT,
        solar_model::{Earth, Moon, OrbitModel, Solarobj, Sun},
    },
    input::{EnvInitData, RuntimeParameters},
    output,
    types::{self, Array3d, LLH},
};

// Hill Sphere source
// Title: Gravitational Spheres of the Major Planets, Moon and Sun
// Journal: Soviet Astronomy, Vol. 7, p.618
// Authors: Chebotarev, G. A.
const AU_METER: f64 = 1.496e+11;
const EARTH_HILL_SPHERE_RADIUS: f64 = 0.01001 * AU_METER;
const LUNAR_HILL_SPHERE_RADIUS: f64 = 58050000.0;

#[derive(Clone)]
pub struct Environment {
    pub start_time: chrono::DateTime<Utc>, // Start time of the simulation as a constant.
    pub sim_time_s: f64,                   // Simulation time in seconds.
    pub step_count: u64,                   // Number of simulation steps run.
    pub current_time: chrono::DateTime<Utc>, // Current time of simulation in UTC.
    // All solar bodies below are synced to future time
    pub sun: Sun,
    pub earth: Earth,
    pub moon: Moon,
}

impl Environment {
    /// Updates the solar system objects within the environment.
    ///
    /// ### Argument
    /// * 'datetime' - New datetime of simulation.
    ///
    fn update_solar_objs(&mut self, datetime: &DateTime<chrono::Utc>) {
        let et = spice_ephemeris_time(datetime);

        self.sun.model.update_position_velocity(et);
        self.earth.model.update_position_velocity(et);
        self.moon.model.update_position_velocity(et);
    }

    /// Main simulation environment advancement method. This method should
    /// perform all required operations to move the simulation environment from
    /// step_count to step_count + 1.
    ///
    /// ### Argument
    /// 'runtime_params" - Simulation parameters gathered at invocation time
    ///
    pub fn advance_simulation_environment(&mut self, runtime_params: &RuntimeParameters) {
        // Advance the simulation environment by configured simulation time step
        self.sim_time_s += runtime_params.sim_time_step as f64;
        self.step_count += 1;
        self.current_time =
            self.start_time + Duration::milliseconds((self.sim_time_s * 1000.0) as i64);

        self.update_solar_objs(&self.current_time.into());

        // Update body model state for use on next sim step.
        self.earth.set_query_space_weather_index(self.current_time);
    }

    // Calculate solar ecliptic coordinates for a given simulation object from
    // the soi inertial reference frame.
    pub fn calculate_helio_coords(&self, sim_obj: &SimobjT) -> Array3d {
        match sim_obj.state.soi {
            Solarobj::Sun => self.sun.model.state.coords + sim_obj.state.coords,
            Solarobj::Earth => {
                self.earth.model.state.coords
                    + bodies::common::equatorial_to_ecliptic(
                        &sim_obj.state.coords,
                        self.earth.attr.obliquity,
                    )
            }
            Solarobj::Moon => {
                self.moon.model.state.coords
                    + bodies::common::equatorial_to_ecliptic(
                        &sim_obj.state.coords,
                        self.moon.attr.obliquity,
                    )
            }
        }
    }

    pub fn calculate_fixed_coords(&self, sim_obj: &SimobjT) -> LLH {
        match sim_obj.state.soi {
            Solarobj::Sun => LLH::default(),
            Solarobj::Earth => self
                .earth
                .eci2geo(&sim_obj.state.coords, days_since_j2000(&self.current_time)),
            Solarobj::Moon => LLH::default(),
        }
    }

    pub fn new(runtime_params: &RuntimeParameters, init_data: EnvInitData) -> Environment {
        let et = spice_ephemeris_time(&runtime_params.date);

        let sun_precalc = Sun::new(et);
        let earth_precalc = Earth::new(et, &init_data.earth_sw);
        let moon_precalc = Moon::new(et);

        let mut env = Environment {
            start_time: runtime_params.date,
            sim_time_s: 0f64,
            step_count: 0u64,
            current_time: runtime_params.date,
            sun: sun_precalc,
            earth: earth_precalc,
            moon: moon_precalc,
        };

        // This forces a hard update which calculates future locations of solar bodies
        env.update_solar_objs(&runtime_params.date);
        env.earth.set_query_space_weather_index(runtime_params.date);

        env
    }

    fn check_is_within_earth_hill_sphere(&self, sim_obj: &mut SimobjT) -> bool {
        let distance_to_earth =
            types::l2_norm(&(sim_obj.state.coord_helio - self.earth.model.state.coords));

        distance_to_earth <= EARTH_HILL_SPHERE_RADIUS
    }

    fn check_is_within_lunar_hill_sphere(&self, sim_obj: &mut SimobjT) -> bool {
        let distance_to_moon =
            types::l2_norm(&(sim_obj.state.coord_helio - self.moon.model.state.coords));

        distance_to_moon <= LUNAR_HILL_SPHERE_RADIUS
    }

    /// If the simulation object is within the hill sphere radius of the solar object in question, switch to that SOI.
    ///
    /// ### Argument
    /// * 'sim_object' - The simulation object to be checked
    ///
    pub fn check_switch_soi(&self, sim_obj: &mut SimobjT) {
        match sim_obj.state.soi {
            // SWitch from Solar SOI into Earth SOI.
            // TODO: Calculations below have not been verified as accurate.
            Solarobj::Sun => {
                if self.check_is_within_earth_hill_sphere(sim_obj) {
                    sim_obj.state.coords = bodies::common::ecliptic_to_equatorial(
                        &(sim_obj.state.coord_helio - self.earth.model.state.coords),
                        self.earth.attr.obliquity,
                    );
                    sim_obj.state.velocity = bodies::common::ecliptic_to_equatorial(
                        &(sim_obj.state.velocity - self.earth.model.state.velocity),
                        self.earth.attr.obliquity,
                    );
                    sim_obj.state.soi = Solarobj::Earth;
                    // Recursive call to handle bodies within switched soi
                    self.check_switch_soi(sim_obj);
                }
            }
            Solarobj::Earth => {
                let sim_distance_to_earth =
                    types::l2_norm(&(sim_obj.state.coord_helio - self.earth.model.state.coords));

                // Switch from Earth SOI into the Solar SOI.
                // Calculations below are verified as accurate.
                if sim_distance_to_earth > EARTH_HILL_SPHERE_RADIUS {
                    sim_obj.state.coords = sim_obj.state.coord_helio - self.sun.model.state.coords;
                    sim_obj.state.velocity = bodies::common::equatorial_to_ecliptic(
                        &sim_obj.state.velocity,
                        self.earth.attr.obliquity,
                    ) + self.earth.model.state.velocity;
                    sim_obj.state.soi = Solarobj::Sun;
                    self.check_switch_soi(sim_obj);
                // Switch from Earth SOI into the Lunar SOI.
                // TODO: Calculations below have not been verified as accurate.
                } else if self.check_is_within_lunar_hill_sphere(sim_obj) {
                    sim_obj.state.coords = bodies::common::ecliptic_to_equatorial(
                        &(sim_obj.state.coord_helio - self.moon.model.state.coords),
                        self.moon.attr.obliquity,
                    );
                    let helio_velocity = bodies::common::equatorial_to_ecliptic(
                        &sim_obj.state.velocity,
                        self.earth.attr.obliquity,
                    ) + self.earth.model.state.velocity;
                    sim_obj.state.velocity = bodies::common::ecliptic_to_equatorial(
                        &(helio_velocity - self.moon.model.state.velocity),
                        self.moon.attr.obliquity,
                    );
                    sim_obj.state.soi = Solarobj::Moon;
                }
            }
            Solarobj::Moon => {
                let sim_distance_to_moon =
                    types::l2_norm(&(sim_obj.state.coord_helio - self.moon.model.state.coords));

                // Switch from Lunar SOI into Earth SOI.
                // TODO: Calculations below have not been verified as accurate.
                if sim_distance_to_moon > LUNAR_HILL_SPHERE_RADIUS {
                    sim_obj.state.coords = bodies::common::ecliptic_to_equatorial(
                        &(sim_obj.state.coord_helio - self.earth.model.state.coords),
                        self.earth.attr.obliquity,
                    );
                    let helio_velocity = bodies::common::equatorial_to_ecliptic(
                        &sim_obj.state.velocity,
                        self.moon.attr.obliquity,
                    ) + self.moon.model.state.velocity;
                    sim_obj.state.velocity = bodies::common::ecliptic_to_equatorial(
                        &(helio_velocity - self.earth.model.state.velocity),
                        self.earth.attr.obliquity,
                    );
                    self.check_switch_soi(sim_obj);
                }
                // Moon has no solar bodies orbiting it
            }
        }
    }

    /// Convert solar data into output form
    ///
    /// ### Return
    /// Returns solar data in form which is used for output
    ///
    pub fn sun_to_output_form(&self) -> output::SolarObjectOut {
        output::SolarObjectOut {
            name: Solarobj::Sun.to_string(),
            sim_time: self.sim_time_s,
            x_coord: self.sun.model.state.coords.x as f32,
            y_coord: self.sun.model.state.coords.y as f32,
            z_coord: self.sun.model.state.coords.z as f32,
            x_velocity: self.sun.model.state.velocity.x,
            y_velocity: self.sun.model.state.velocity.y,
            z_velocity: self.sun.model.state.velocity.z,
        }
    }

    /// Convert earth data into output form
    ///
    /// ### Return
    /// Returns earth data in form which is used for output
    ///
    pub fn earth_to_output_form(&self) -> output::SolarObjectOut {
        output::SolarObjectOut {
            name: Solarobj::Earth.to_string(),
            sim_time: self.sim_time_s,
            x_coord: self.earth.model.state.coords.x as f32,
            y_coord: self.earth.model.state.coords.y as f32,
            z_coord: self.earth.model.state.coords.z as f32,
            x_velocity: self.earth.model.state.velocity.x,
            y_velocity: self.earth.model.state.velocity.y,
            z_velocity: self.earth.model.state.velocity.z,
        }
    }

    /// Convert lunar data into output form
    ///
    /// ### Return
    /// Return lunar data in form which is used for output
    ///
    pub fn moon_to_output_form(&self) -> output::SolarObjectOut {
        output::SolarObjectOut {
            name: Solarobj::Moon.to_string(),
            sim_time: self.sim_time_s,
            x_coord: self.moon.model.state.coords.x as f32,
            y_coord: self.moon.model.state.coords.y as f32,
            z_coord: self.moon.model.state.coords.z as f32,
            x_velocity: self.moon.model.state.velocity.x,
            y_velocity: self.moon.model.state.velocity.y,
            z_velocity: self.moon.model.state.velocity.z,
        }
    }
}
