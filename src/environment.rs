use chrono::{DateTime, Duration, Utc};

use crate::{
    bodies::{
        self,
        common::days_since_j2000,
        sim_object::SimobjT,
        solar_model::{make_earth, make_moon, make_sun, Earth, KeplerModel, Moon, Solarobj, Sun},
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

pub struct Environment {
    pub start_time: chrono::DateTime<Utc>, // Start time of the simulation as a constant.
    sim_time_s: f64,                       // Simulation time in seconds.
    pub current_time: chrono::DateTime<Utc>, // Current time of simulation in UTC.
    last_day_update_s: f64,                // Last time a hard simulation model update was done.
    future_day_update_s: f64, // Next time a hard simulation model update should be done.
    // All solar bodies below are synced to future time
    pub sun: Sun,
    pub earth: Earth,
    pub moon: Moon,
}

impl Environment {
    /// Updates the solar system objects within the environment.
    ///
    /// ### Argument
    /// * 'up_day' - New datetime of simulation.
    ///
    fn update_solar_objs(&mut self, up_day: &DateTime<chrono::Utc>) {
        let new_day = days_since_j2000(up_day);

        // Update each solar body within simulation
        self.sun.model.state.coords.ahead_coords =
            self.sun.model.ecliptic_cartesian_coords(new_day);
        self.earth.model.state.coords.ahead_coords =
            self.earth.model.ecliptic_cartesian_coords(new_day);
        // Calculate new location for moon and convert to heliocentric coords
        self.moon.model.state.coords.ahead_coords =
            self.moon.model.ecliptic_cartesian_coords(new_day)
                + self.earth.model.state.coords.ahead_coords;
    }

    /// Hard update on all solar objects within the simulation.
    /// Simulation time fields within environment should be set accordingly
    /// before this method is called.
    ///
    /// ### Argument
    /// 'runtime_params" - Simulation parameters gathered at invocation time
    ///
    fn update(&mut self, runtime_params: &RuntimeParameters) {
        self.last_day_update_s = self.sim_time_s;

        self.update_solar_objs(&self.current_time.into());

        // Set the corresponding last and current coords to ahead time as ahead time is now current
        self.sun.model.state.coords.behind_coords = self.sun.model.state.coords.ahead_coords;
        self.sun.model.state.coords.current_coords = self.sun.model.state.coords.ahead_coords;

        self.earth.model.state.coords.behind_coords = self.earth.model.state.coords.ahead_coords;
        self.earth.model.state.coords.current_coords = self.earth.model.state.coords.ahead_coords;

        self.moon.model.state.coords.behind_coords = self.moon.model.state.coords.ahead_coords;
        self.moon.model.state.coords.current_coords = self.moon.model.state.coords.ahead_coords;

        self.future_day_update_s = self.sim_time_s + (runtime_params.sim_solar_step as f64);

        // Calculate new positions at future time
        let future_time = self.start_time + Duration::seconds(self.future_day_update_s as i64);
        self.update_solar_objs(&future_time);

        self.sun.model.state.velocity = self
            .sun
            .model
            .state
            .coords
            .calc_velocity_full_range(runtime_params.sim_solar_step as f64);
        self.earth.model.state.velocity = self
            .earth
            .model
            .state
            .coords
            .calc_velocity_relative_full_range(
                &self.sun.model.state.coords,
                runtime_params.sim_solar_step as f64,
            );
        self.moon.model.state.velocity = self
            .moon
            .model
            .state
            .coords
            .calc_velocity_relative_full_range(
                &self.earth.model.state.coords,
                runtime_params.sim_solar_step as f64,
            );
    }

    /// Advance the simulation by the simulation step time. This function will force a hard update on
    /// all solar objects within the simulation if current simulation time exceeds or is equal to future
    /// simulation time. In between hard environment updates, positions of all solar objects are linearly interpolated.
    ///
    /// ### Argument
    /// 'runtime_params" - Simulation parameters gathered at invocation time
    ///
    pub fn advance_simulation_environment(&mut self, runtime_params: &RuntimeParameters) {
        // Advance the simulation environment by configured simulation time step
        self.sim_time_s += runtime_params.sim_time_step as f64;
        self.current_time =
            self.start_time + Duration::milliseconds((self.sim_time_s * 1000.0) as i64);
        // Perform a hard update if advanced current simulation time would exceed future simulation time
        if self.sim_time_s >= self.future_day_update_s {
            self.update(runtime_params);
            // Exit out of function as hard update handled advance of simulation environment
            return;
        }

        // Perform a linear interpolation from the current solar object coords to the future solar coords
        let interp_point = (self.sim_time_s - self.last_day_update_s)
            / (self.future_day_update_s - self.last_day_update_s);

        // Linear interpolate for all solar objs within simulation
        self.sun.model.state.coords.lerp_set(interp_point);
        self.earth.model.state.coords.lerp_set(interp_point);
        self.moon.model.state.coords.lerp_set(interp_point);
    }

    // Calculate solar ecliptic coordinates for a given simulation object from
    // the soi inertial reference frame.
    pub fn calculate_se_coords(&self, sim_obj: &SimobjT) -> Array3d {
        match sim_obj.soi {
            Solarobj::Sun => self.sun.model.state.coords.current_coords + sim_obj.coords,
            Solarobj::Earth => {
                self.earth.model.state.coords.current_coords
                    + bodies::common::equatorial_to_ecliptic(
                        &sim_obj.coords,
                        self.earth.attr.obliquity,
                    )
            }
            Solarobj::Moon => {
                self.moon.model.state.coords.current_coords
                    + bodies::common::equatorial_to_ecliptic(
                        &sim_obj.coords,
                        self.moon.attr.obliquity,
                    )
            }
        }
    }

    pub fn calculate_fixed_coords(&self, sim_obj: &SimobjT) -> LLH {
        match sim_obj.soi {
            Solarobj::Sun => LLH::default(),
            Solarobj::Earth => self
                .earth
                .eci2geo(&sim_obj.coords, days_since_j2000(&self.current_time)),
            Solarobj::Moon => LLH::default(),
        }
    }

    pub fn new(runtime_params: &RuntimeParameters, init_data: EnvInitData) -> Environment {
        let day = days_since_j2000(&runtime_params.date);

        let sun_precalc = make_sun();
        let earth_precalc = make_earth(day, &init_data.earth_sw);
        let moon_precalc = make_moon(day, &earth_precalc.model.state.coords.current_coords);

        let mut env = Environment {
            start_time: runtime_params.date,
            sim_time_s: 0f64,
            current_time: runtime_params.date,
            last_day_update_s: 0f64,
            future_day_update_s: runtime_params.sim_solar_step as f64,
            sun: sun_precalc,
            earth: earth_precalc,
            moon: moon_precalc,
        };

        // This forces a hard update which calculates future locations of solar bodies
        env.update(runtime_params);

        env
    }

    /// Get simulation time in seconds
    pub fn get_sim_time(&self) -> f64 {
        self.sim_time_s
    }

    fn check_is_within_earth_hill_sphere(&self, sim_obj: &mut SimobjT) -> bool {
        let distance_to_earth =
            types::l2_norm(&(sim_obj.coords_abs - self.earth.model.state.coords.current_coords));

        distance_to_earth <= EARTH_HILL_SPHERE_RADIUS
    }

    fn check_is_within_lunar_hill_sphere(&self, sim_obj: &mut SimobjT) -> bool {
        let distance_to_moon =
            types::l2_norm(&(sim_obj.coords_abs - self.moon.model.state.coords.current_coords));

        distance_to_moon <= LUNAR_HILL_SPHERE_RADIUS
    }

    fn switch_in_soi(
        sim_obj: &mut SimobjT,
        to_soi: Solarobj,
        soi_coords: &Array3d,
        soi_velocity: &Array3d,
    ) {
        sim_obj.coords = sim_obj.coords - soi_coords;
        sim_obj.velocity = sim_obj.velocity - soi_velocity;
        sim_obj.soi = to_soi;
    }

    fn switch_out_soi(
        sim_obj: &mut SimobjT,
        to_soi: Solarobj,
        soi_coords: &Array3d,
        soi_velocity: &Array3d,
    ) {
        sim_obj.coords = sim_obj.coords + soi_coords;
        sim_obj.velocity = sim_obj.velocity + soi_velocity;
        sim_obj.soi = to_soi;
    }

    /// If the simulation object is within the hill sphere radius of the solar object in question, switch to that SOI.
    ///
    /// ### Argument
    /// * 'sim_object' - The simulation object to be checked
    ///
    pub fn check_switch_soi(&self, sim_obj: &mut SimobjT) {
        match sim_obj.soi {
            Solarobj::Sun => {
                if self.check_is_within_earth_hill_sphere(sim_obj) {
                    Self::switch_in_soi(
                        sim_obj,
                        Solarobj::Earth,
                        &self.earth.model.state.coords.current_coords,
                        &self.earth.model.state.velocity,
                    );
                    // Recursive call to handle bodies within switched soi
                    self.check_switch_soi(sim_obj);
                }
            }
            Solarobj::Earth => {
                let sim_distance_to_earth = types::l2_norm(
                    &(sim_obj.coords_abs - self.earth.model.state.coords.current_coords),
                );

                // If the simulation object is outside of the earth hill sphere set the soi to solar
                if sim_distance_to_earth > EARTH_HILL_SPHERE_RADIUS {
                    Self::switch_out_soi(
                        sim_obj,
                        Solarobj::Sun,
                        &(self.earth.model.state.coords.current_coords
                            - self.sun.model.state.coords.current_coords),
                        &self.earth.model.state.velocity,
                    );
                    // Recursive call to handle bodies within switched soi, used for overshoots/teleportation
                    self.check_switch_soi(sim_obj);
                } else if self.check_is_within_lunar_hill_sphere(sim_obj) {
                    Self::switch_in_soi(
                        sim_obj,
                        Solarobj::Moon,
                        &self.moon.model.state.coords.current_coords,
                        &self.moon.model.state.velocity,
                    )
                }
            }
            Solarobj::Moon => {
                let sim_distance_to_moon = types::l2_norm(
                    &(sim_obj.coords_abs - self.moon.model.state.coords.current_coords),
                );

                // If the simulation object is outside the lunar hill sphere set the soi to earth
                if sim_distance_to_moon > LUNAR_HILL_SPHERE_RADIUS {
                    Self::switch_out_soi(
                        sim_obj,
                        Solarobj::Earth,
                        &(self.moon.model.state.coords.current_coords
                            - self.earth.model.state.coords.current_coords),
                        &self.moon.model.state.velocity,
                    );
                    // Recursive call to handle bodies within switched soi, used for overshoots/teleportation
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
            x_coord: self.sun.model.state.coords.current_coords.x as f32,
            y_coord: self.sun.model.state.coords.current_coords.y as f32,
            z_coord: self.sun.model.state.coords.current_coords.z as f32,
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
            x_coord: self.earth.model.state.coords.current_coords.x as f32,
            y_coord: self.earth.model.state.coords.current_coords.y as f32,
            z_coord: self.earth.model.state.coords.current_coords.z as f32,
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
            x_coord: self.moon.model.state.coords.current_coords.x as f32,
            y_coord: self.moon.model.state.coords.current_coords.y as f32,
            z_coord: self.moon.model.state.coords.current_coords.z as f32,
            x_velocity: self.moon.model.state.velocity.x,
            y_velocity: self.moon.model.state.velocity.y,
            z_velocity: self.moon.model.state.velocity.z,
        }
    }
}
