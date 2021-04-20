use crate::{input::SimulationParameters, output, types};
use chrono::{DateTime, Duration, TimeZone, Utc};
use serde::{Deserialize, Serialize};
use strum_macros::Display;
use types::Array3d;

const METERS_PER_ASTRONOMICAL_UNIT: f64 = 1.4959787e+11;
const METERS_PER_EARTH_EQUATORIAL_RADIUS: f64 = 6378140.0;
const EARTH_RADII_PER_ASTRONOMICAL_UNIT: f64 =
    METERS_PER_ASTRONOMICAL_UNIT / METERS_PER_EARTH_EQUATORIAL_RADIUS; // 23454.78
const AU_METER: f64 = 1.496e+11;

// Hill Sphere source
// Title: Gravitational Spheres of the Major Planets, Moon and Sun
// Journal: Soviet Astronomy, Vol. 7, p.618
// Authors: Chebotarev, G. A.
const EARTH_HILL_SPHERE_RADIUS: f64 = 0.01001 * AU_METER;
const LUNAR_HILL_SPHERE_RADIUS: f64 = 58050000.0;

#[derive(Serialize, Deserialize)]
pub struct InitData {
    pub date: String,             // Datetime in ISO 8601 format
    pub debris: Vec<SimobjT>,     // Debris objects
    pub spacecraft: Vec<SimobjT>, // Spacecraft objects
}

#[derive(Serialize)]
pub enum SimObjectType {
    Spacecraft,
    Debris,
}

impl Default for SimObjectType {
    fn default() -> Self {
        SimObjectType::Spacecraft
    }
}

#[derive(Serialize, Deserialize)]
pub struct SimobjT {
    #[serde(skip_deserializing)]
    pub id: u32,
    #[serde(skip_deserializing)]
    pub sim_object_type: SimObjectType,
    #[serde(skip_deserializing)]
    pub coords_abs: Array3d,
    pub soi: Solarobj,
    pub coords: Array3d,
    pub velocity: Array3d,
    pub drag_area: f64,
    pub mass: f64,
}

impl SimobjT {
    fn type_of(&self) -> String {
        match self.sim_object_type {
            SimObjectType::Spacecraft => String::from("Spacecraft"),
            SimObjectType::Debris => String::from("Debris"),
        }
    }

    pub fn to_output_form(&self, sim_time: f64) -> output::SimulationObjectParameters {
        let abs_coords = &self.coords_abs;
        let coords = &self.coords;
        let velocity = &self.velocity;

        output::SimulationObjectParameters {
            id: self.id,
            sim_time,
            soi: self.soi.to_string(),
            x_abs_coord: abs_coords.x,
            y_abs_coord: abs_coords.y,
            z_abs_coord: abs_coords.z,
            x_coord: coords.x,
            y_coord: coords.y,
            z_coord: coords.z,
            x_velocity: velocity.x,
            y_velocity: velocity.y,
            z_velocity: velocity.z,
        }
    }
}

pub struct Environment {
    last_day_update_s: f64,
    sim_time_s: f64, // Simulation time in seconds
    future_day_update_s: f64,
    pub start_time: chrono::DateTime<Utc>,
    // All solar bodies below are synced to future time
    pub sun: Sun,
    pub earth: Earth,
    pub moon: PlanetPS,
}

impl Environment {
    /// Updates the solar system objects within the environment.
    ///
    /// ### Argument
    /// * 'up_day' - New datetime of simulation.
    ///
    fn update_solar_objs(&mut self, up_day: &DateTime<chrono::Utc>) {
        let new_day = Self::datetime_to_days(up_day);

        // Update each solar body within simulation
        self.sun.coords.ahead_coords = self.sun.ecliptic_cartesian_coords(new_day);
        self.earth.coords.ahead_coords = self.earth.ecliptic_cartesian_coords(new_day);
        // Calculate new location for moon and convert to heliocentric coords
        self.moon.coords.ahead_coords =
            self.moon.ecliptic_cartesian_coords(new_day) + self.earth.coords.ahead_coords;
    }

    /// Hard update on all solar objects within the simulation.
    ///
    /// ### Argument
    /// 'sim_params" - Simulation parameters gathered at invocation time
    ///
    fn update(&mut self, sim_params: &SimulationParameters) {
        self.last_day_update_s = self.sim_time_s;

        // Calculate postions at current time (postions at future time, TODO optimise)
        let new_time = self.start_time + Duration::seconds(self.sim_time_s as i64);
        self.update_solar_objs(&new_time);

        // Set the corresponding last and current coords to ahead time as ahead time is now current
        self.sun.coords.behind_coords = self.sun.coords.ahead_coords;
        self.sun.coords.current_coords = self.sun.coords.ahead_coords;

        self.earth.coords.behind_coords = self.earth.coords.ahead_coords;
        self.earth.coords.current_coords = self.earth.coords.ahead_coords;

        self.moon.coords.behind_coords = self.moon.coords.ahead_coords;
        self.moon.coords.current_coords = self.moon.coords.ahead_coords;

        self.future_day_update_s = self.sim_time_s + (sim_params.sim_solar_step as f64);

        // Calculate new positions at future time
        let future_time = self.start_time + Duration::seconds(self.future_day_update_s as i64);
        self.update_solar_objs(&future_time);

        self.sun.velocity = self
            .sun
            .coords
            .calc_velocity_full_range(sim_params.sim_solar_step as f64);
        self.earth.velocity = self
            .earth
            .coords
            .calc_velocity_relative_full_range(&self.sun.coords, sim_params.sim_solar_step as f64);
        self.moon.velocity = self.moon.coords.calc_velocity_relative_full_range(
            &self.earth.coords,
            sim_params.sim_solar_step as f64,
        );
    }

    /// Advance the simulation by the simulation step time. This function will force a hard update on
    /// all solar objects within the simulation if current simulation time exceeds or is equal to future
    /// simulation time. In between hard environment updates, positions of all solar objects are linearly interpolated.
    ///
    /// ### Argument
    /// 'sim_params" - Simulation parameters gathered at invocation time
    ///
    pub fn advance_simulation_environment(&mut self, sim_params: &SimulationParameters) {
        // Advance the simulation environment by configured simulation time step
        self.sim_time_s += sim_params.sim_time_step as f64;
        // Perform a hard update if advanced current simulation time would exceed future simulation time
        if self.sim_time_s >= self.future_day_update_s {
            self.update(sim_params);
            // Exit out of function as hard update handled advance of simulation environment
            return;
        }

        // Perform a linear interpolation from the current solar object coords to the future solar coords
        let interp_point = (self.sim_time_s - self.last_day_update_s)
            / (self.future_day_update_s - self.last_day_update_s);

        // Linear interpolate for all solar objs within simulation
        self.sun.coords.lerp_set(interp_point);
        self.earth.coords.lerp_set(interp_point);
        self.moon.coords.lerp_set(interp_point);
    }

    /// Calculate the absolute coordinates for the simulation object
    ///
    /// ### Arguments
    /// * 'sim_obj' - The simulation object
    ///
    /// ### Return
    /// The absolute coordinates for the input simulation object.
    ///
    pub fn calculate_abs_coords(&self, sim_obj: &SimobjT) -> Array3d {
        match sim_obj.soi {
            Solarobj::Sun { attr: _, .. } => self.sun.coords.current_coords + sim_obj.coords,
            Solarobj::Earth { attr: _, .. } => self.earth.coords.current_coords + sim_obj.coords,
            // Lunar coords are abs also
            Solarobj::Moon { attr: _, .. } => self.moon.coords.current_coords + sim_obj.coords,
        }
    }

    /// Calculates a delta for provided datetime from 0/Jan/2000 00:00 UTC
    ///
    /// ### Argument
    /// * 'datetime' - User provided datetime object.
    ///
    /// ### Return
    ///     The delta from 30/Dec/1999 20:00:00 UTC in days.
    ///
    fn datetime_to_days(datetime_obj: &DateTime<chrono::Utc>) -> f64 {
        let origin_dt = chrono::Utc.ymd(1999, 12, 30).and_hms(20, 0, 0);

        (*datetime_obj - origin_dt).num_seconds() as f64 / 86400f64
    }

    /// Creates the initial solar system
    ///
    /// ### Argument
    /// * 'day' - The day value greater than zero. From 2000-01-01
    ///
    /// ### Return
    ///     new Environment loaded with all possible solar system objects tracked by the
    ///     simulation.
    ///
    pub fn new(start_time: DateTime<Utc>, sim_params: &SimulationParameters) -> Environment {
        let day = Environment::datetime_to_days(&start_time);

        let sun_precalc = make_sun();
        let earth_precalc = make_earth(day);
        let moon_precalc = make_moon(day, &earth_precalc.coords.current_coords);

        let mut env = Environment {
            last_day_update_s: 0f64,
            sim_time_s: 0f64,
            future_day_update_s: sim_params.sim_solar_step as f64,
            start_time,
            sun: sun_precalc,
            earth: earth_precalc,
            moon: moon_precalc,
        };

        // This forces a hard update which calculates future locations of solar bodies
        env.update(sim_params);

        env
    }

    /// Get simulation time in seconds
    pub fn get_sim_time(&self) -> f64 {
        self.sim_time_s
    }

    fn check_is_within_earth_hill_sphere(&self, sim_obj: &mut SimobjT) -> bool {
        let distance_to_earth =
            types::l2_norm(&(sim_obj.coords_abs - self.earth.coords.current_coords));

        distance_to_earth <= EARTH_HILL_SPHERE_RADIUS
    }

    fn check_is_within_lunar_hill_sphere(&self, sim_obj: &mut SimobjT) -> bool {
        let distance_to_moon =
            types::l2_norm(&(sim_obj.coords_abs - self.moon.coords.current_coords));

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
            Solarobj::Sun { attr: _ } => {
                if self.check_is_within_earth_hill_sphere(sim_obj) {
                    Self::switch_in_soi(
                        sim_obj,
                        Solarobj::Earth { attr: None },
                        &self.earth.coords.current_coords,
                        &self.earth.velocity,
                    );
                    // Recursive call to handle bodies within switched soi
                    self.check_switch_soi(sim_obj);
                }
            }
            Solarobj::Earth { attr: _ } => {
                let sim_distance_to_earth =
                    types::l2_norm(&(sim_obj.coords_abs - self.earth.coords.current_coords));

                // If the simulation object is outside of the earth hill sphere set the soi to solar
                if sim_distance_to_earth > EARTH_HILL_SPHERE_RADIUS {
                    Self::switch_out_soi(
                        sim_obj,
                        Solarobj::Sun { attr: None },
                        &(self.earth.coords.current_coords - self.sun.coords.current_coords),
                        &self.earth.velocity,
                    );
                    // Recursive call to handle bodies within switched soi, used for overshoots/teleportation
                    self.check_switch_soi(sim_obj);
                } else if self.check_is_within_lunar_hill_sphere(sim_obj) {
                    Self::switch_in_soi(
                        sim_obj,
                        Solarobj::Moon { attr: None },
                        &self.moon.coords.current_coords,
                        &self.moon.velocity,
                    )
                }
            }
            Solarobj::Moon { attr: _ } => {
                let sim_distance_to_moon =
                    types::l2_norm(&(sim_obj.coords_abs - self.moon.coords.current_coords));

                // If the simulation object is outside the lunar hill sphere set the soi to earth
                if sim_distance_to_moon > LUNAR_HILL_SPHERE_RADIUS {
                    Self::switch_out_soi(
                        sim_obj,
                        Solarobj::Earth { attr: None },
                        &(self.moon.coords.current_coords - self.earth.coords.current_coords),
                        &self.moon.velocity,
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
            name: self.sun.get_solar_object().to_string(),
            sim_time: self.sim_time_s,
            x_coord: self.sun.coords.current_coords.x as f32,
            y_coord: self.sun.coords.current_coords.y as f32,
            z_coord: self.sun.coords.current_coords.z as f32,
            x_velocity: self.sun.velocity.x,
            y_velocity: self.sun.velocity.y,
            z_velocity: self.sun.velocity.z,
        }
    }

    /// Convert earth data into output form
    ///
    /// ### Return
    /// Returns earth data in form which is used for output
    ///
    pub fn earth_to_output_form(&self) -> output::SolarObjectOut {
        output::SolarObjectOut {
            name: self.earth.get_solar_object().to_string(),
            sim_time: self.sim_time_s,
            x_coord: self.earth.coords.current_coords.x as f32,
            y_coord: self.earth.coords.current_coords.y as f32,
            z_coord: self.earth.coords.current_coords.z as f32,
            x_velocity: self.earth.velocity.x,
            y_velocity: self.earth.velocity.y,
            z_velocity: self.earth.velocity.z,
        }
    }

    /// Convert lunar data into output form
    ///
    /// ### Return
    /// Return lunar data in form which is used for output
    ///
    pub fn moon_to_output_form(&self) -> output::SolarObjectOut {
        output::SolarObjectOut {
            name: self.moon.get_solar_object().to_string(),
            sim_time: self.sim_time_s,
            x_coord: self.moon.coords.current_coords.x as f32,
            y_coord: self.moon.coords.current_coords.y as f32,
            z_coord: self.moon.coords.current_coords.z as f32,
            x_velocity: self.moon.velocity.x,
            y_velocity: self.moon.velocity.y,
            z_velocity: self.moon.velocity.z,
        }
    }
}

#[derive(Display, Clone, Deserialize, Serialize)]
#[serde(tag = "type")]
pub enum Solarobj {
    #[strum(serialize = "sun")]
    Sun { attr: Option<SolarAttr> },
    #[strum(serialize = "earth")]
    Earth { attr: Option<SolarAttr> },
    #[strum(serialize = "moon")]
    Moon { attr: Option<SolarAttr> },
}

#[derive(Clone, Deserialize, Serialize)]
pub struct SolarAttr {
    radius: f64, // meters
    mass: f64,   // kg
}

impl Solarobj {
    pub fn get_mass_kg(&self) -> f64 {
        const ERROR_MSG: &str = "Enum has some field for attr";

        match self {
            Solarobj::Sun { attr } => attr.as_ref().expect(ERROR_MSG).mass,
            Solarobj::Earth { attr } => attr.as_ref().expect(ERROR_MSG).mass,
            Solarobj::Moon { attr } => attr.as_ref().expect(ERROR_MSG).mass,
        }
    }
}

#[derive(Clone)]
pub struct PlanetPS {
    // See  http://www.stjarnhimlen.se/comp/ppcomp.html#4
    solartype: Solarobj, // Type enum of the solar obj
    pub coords: SolarobjCoords,
    pub velocity: Array3d,
    n0: f64,
    nc: f64, // N0 = longitude of the ascending node (deg).  Nc = rate of change in deg/day
    i0: f64,
    ic: f64, // inclination to the ecliptic (deg)
    w0: f64,
    wc: f64, // argument of perihelion (deg)
    a0: f64,
    ac: f64, // semi-major axis, or mean distance from Sun (AU)
    e0: f64,
    ec: f64, // eccentricity (0=circle, 0..1=ellipse, 1=parabola)
    m0: f64,
    mc: f64, // M0 = mean anomaly  (deg) (0 at perihelion; increases uniformly with time).  Mc ("mean motion") = rate of change
    mag_base: f64,
    mag_phase_factor: f64,
    mag_nonlinear_factor: f64,
    mag_nonlinear_exponent: f64,
}

#[derive(Clone)]
pub struct Earth {
    solartype: Solarobj,
    pub coords: SolarobjCoords,
    pub velocity: Array3d,
}

#[derive(Clone)]
pub struct Sun {
    solartype: Solarobj,
    pub coords: SolarobjCoords,
    pub velocity: Array3d,
}

/// Provides utilities for calculating planetary bodies with a Kepler model
mod kepler_utilities {
    use crate::{
        bodies::{PlanetPS, EARTH_RADII_PER_ASTRONOMICAL_UNIT},
        types::Array3d,
    };
    use std::f64::{self, consts};

    /// Calculate the eccentric anomaly for a given body.
    /// ### Arguments
    /// * 'e' - TODO
    /// * 'm' - TODO
    ///
    /// ### Returns
    ///      The eccentric anomaly for the provided input parameters.
    pub fn eccentric_anomaly(e: f64, m: f64) -> f64 {
        let deg_from_rad = 180f64 / consts::PI;
        let mut ecc: f64 = m + (e * sin_deg!(m) * (1f64 + (e * cos_deg!(m))));

        let mut iters = 0;
        loop {
            let f: f64 =
                ecc - (ecc - (deg_from_rad * e * sin_deg!(ecc)) - m) / (1f64 - e * cos_deg!(ecc));
            let error = (f - ecc).abs();
            ecc = f;

            if error < 0.1e-8 || iters > 30 {
                break;
            }

            iters += 1;
        }

        ecc
    }

    /// Calculates the mean anomaly for the Sun.
    fn mean_anomaly_of_sun(day: f64) -> f64 {
        (356.0470 + (0.9856002585 * day)) as f64
    }

    /// Calculates the argument of perihelion for the Sun.
    fn sun_argument_of_perihelion(day: f64) -> f64 {
        (282.9404 + (4.70935e-5 * day)) as f64
    }

    /// Calculates the ecliptic latitude and longitude for the given inputs.
    ///
    /// ### Arguments
    /// * 'xh' - Cartesian coordinate in x dimension.
    /// * 'yh' - Cartesian coordinate in y dimension.
    /// * 'zh' - Cartesian coordinate in z dimension.
    ///
    /// ### Return
    ///      The latitude and longitude as a tuple.
    fn ecliptic_lat_lon(xh: f64, yh: f64, zh: f64) -> (f64, f64) {
        (
            atan2_deg!(yh, xh),
            atan2_deg!(zh, (xh * xh + yh * yh).sqrt()),
        )
    }

    pub fn lunar_pertub(body: &PlanetPS, xh: f64, yh: f64, zh: f64, day: f64) -> Array3d {
        let ms = mean_anomaly_of_sun(day); // mean anomaly of Sun
        let ws = sun_argument_of_perihelion(day); // Sun's argument of perihelion
        let ls = ms + ws; // mean longitude of Sun

        let mm = body.mean_anomaly(day); // Moon's mean anomaly
        let nm = body.node_longitude(day); // longitude of Moon's node
        let wm = body.perihelion(day); // Moon's argument of perihelion
        let lm = mm + wm + nm; // Mean longitude of the Moon

        let d = lm - ls; // mean elongation of the Moon
        let f = lm - nm; // argument of latitude for the Moon

        let delta_long = -1.274 * sin_deg!(mm - 2f64*d)       +   // the Evection
            0.658 * sin_deg!(2f64*d)            -   // the Variation
            0.186 * sin_deg!(ms)                -   // the Yearly Equation
            0.059 * sin_deg!(2f64*mm - 2f64*d)  -
            0.057 * sin_deg!(mm - 2f64*d + ms)  +
            0.053 * sin_deg!(mm + 2f64*d)       +
            0.046 * sin_deg!(2f64*d - ms)       +
            0.041 * sin_deg!(mm - ms)           -
            0.035 * sin_deg!(d)                 -   // the Parallactic Equation
            0.031 * sin_deg!(mm + ms)           -
            0.015 * sin_deg!(2f64*f - 2f64*d)
            + 0.011 * sin_deg!(mm - 4f64 * d);

        let delta_lat = -0.173 * sin_deg!(f - 2f64 * d)
            - 0.055 * sin_deg!(mm - f - 2f64 * d)
            - 0.046 * sin_deg!(mm + f - 2f64 * d)
            + 0.033 * sin_deg!(f + 2f64 * d)
            + 0.017 * sin_deg!(2f64 * mm + f);

        let delta_radius = -0.58 * cos_deg!(mm - 2f64 * d) - 0.46 * cos_deg!(2f64 * d);

        let (mut lonecl, mut latecl) = ecliptic_lat_lon(xh, yh, zh);

        let mut r = (xh * xh + yh * yh + zh * zh).sqrt();

        lonecl += delta_long;
        latecl += delta_lat;
        r += delta_radius / EARTH_RADII_PER_ASTRONOMICAL_UNIT;

        let coslon = cos_deg!(lonecl);
        let sinlon = sin_deg!(lonecl);
        let coslat = cos_deg!(latecl);
        let sinlat = sin_deg!(latecl);

        let xp = r * coslon * coslat;
        let yp = r * sinlon * coslat;
        let zp = r * sinlat;

        Array3d {
            x: xp as f64,
            y: yp as f64,
            z: zp as f64,
        }
    }
}

pub trait KeplerModel {
    fn ecliptic_cartesian_coords(&self, day: f64) -> Array3d;

    fn perturb(&self, x: f64, y: f64, z: f64, _day: f64) -> Array3d {
        Array3d { x, y, z }
    }

    fn get_solar_object(&self) -> &Solarobj;
}
#[derive(Clone, Copy)]
pub struct SolarobjCoords {
    pub ahead_coords: Array3d,
    pub current_coords: Array3d, // Initial coords
    pub behind_coords: Array3d,  // Initial coords
}

impl SolarobjCoords {
    /// Calculate linear interpolation of behind and ahead time at interp_point.
    /// Sets the current coords to the result.
    ///
    /// ### Arguments
    /// * 'interp_point' - Point between 0-1 where the interp is set
    ///
    fn lerp_set(&mut self, interp_point: f64) {
        // Precise method, which guarantees v = v1 when t = 1.
        let interp_method =
            |v0: &Array3d, v1: &Array3d, t: f64| -> Array3d { ((1f64 - t) * v0) + (t * v1) };

        self.current_coords = interp_method(&self.behind_coords, &self.ahead_coords, interp_point);
    }

    /// Calculate the velocity of the solar obj from behind to ahead coords.
    ///
    /// ### Arguments
    /// * 'step_time_s' - Sim solar step time in seconds
    ///
    /// ### Returns
    /// The velocity of the solar obj
    ///
    fn calc_velocity_full_range(&self, step_time_s: f64) -> Array3d {
        (self.ahead_coords - self.behind_coords) / step_time_s
    }

    /// Calculate the velocity of the solar obj from behind to ahead coords.
    ///
    /// ### Arguments
    /// * 'step_time_s' - Sim solar step time in seconds
    /// * 'relative' - Solar obj relative coords
    ///
    /// ### Returns
    /// The velocity of the solar obj
    ///
    fn calc_velocity_relative_full_range(
        &self,
        relative: &SolarobjCoords,
        step_time_s: f64,
    ) -> Array3d {
        ((self.ahead_coords - relative.ahead_coords)
            - (self.behind_coords - relative.behind_coords))
            / step_time_s
    }
}

impl Default for SolarobjCoords {
    fn default() -> Self {
        Self {
            ahead_coords: Array3d::default(),
            current_coords: Array3d::default(),
            behind_coords: Array3d::default(),
        }
    }
}

impl KeplerModel for PlanetPS {
    fn ecliptic_cartesian_coords(&self, day: f64) -> Array3d {
        // Default impl
        let a = self.a0 + (day * self.ac);
        let e = self.e0 + (day * self.ec);
        let m_u = self.m0 + (day * self.mc);
        let n_u = self.n0 + (day * self.nc);
        let w = self.w0 + (day * self.wc);
        let i = self.i0 + (day * self.ic);
        let ecc = kepler_utilities::eccentric_anomaly(e, m_u);

        let xv = a * (cos_deg!(ecc) - e);
        let yv = a * ((1.0f64 - e * e).sqrt() * sin_deg!(ecc));

        let v = atan2_deg!(yv, xv); // True anomaly in degrees: the angle from perihelion of the body as seen by the Sun.
        let r = (xv * xv + yv * yv).sqrt(); // Distance from the Sun to the planet in AU

        let cos_n = cos_deg!(n_u);
        let sin_n = sin_deg!(n_u);
        let cosi = cos_deg!(i);
        let sini = sin_deg!(i);
        let cos_vw = cos_deg!(v + w);
        let sin_vw = sin_deg!(v + w);

        // Now we are ready to calculate (unperturbed) ecliptic cartesian heliocentric coordinates.
        let mut xh = r * (cos_n * cos_vw - sin_n * sin_vw * cosi);
        let mut yh = r * (sin_n * cos_vw + cos_n * sin_vw * cosi);
        let mut zh = r * sin_vw * sini;

        // Convert AU to Meters
        xh *= AU_METER;
        yh *= AU_METER;
        zh *= AU_METER;

        self.perturb(xh as f64, yh as f64, zh as f64, day)
    }

    /// Calculates additional perturbations on top of main heliocentric position calculation.
    /// Matches PlanetPS bodies using the type enum.
    ///  
    /// ### Arguments
    ///  * 'x' - X coord
    ///  * 'y' - Y coord
    ///  * 'z' - Z coord
    ///  * 'day' - Day value
    ///
    /// ### Returns
    ///      Cartesian coords with the added perturbations.
    ///
    fn perturb(&self, x: f64, y: f64, z: f64, day: f64) -> Array3d {
        match &self.solartype {
            Solarobj::Moon { attr: _ } => kepler_utilities::lunar_pertub(self, x, y, z, day),
            _ => Array3d { x, y, z },
        }
    }

    fn get_solar_object(&self) -> &Solarobj {
        &self.solartype
    }
}

impl PlanetPS {
    fn mean_anomaly(&self, day: f64) -> f64 {
        self.m0 + (day * self.mc)
    }

    fn node_longitude(&self, day: f64) -> f64 {
        self.n0 + (day * self.nc)
    }

    fn perihelion(&self, day: f64) -> f64 {
        self.w0 + (day * self.wc)
    }
}

impl KeplerModel for Earth {
    /// Calculate the position of Earth relative to the Sun.
    /// Calls function earth_ecliptic_cartesian_coords in kepler_utilities
    ///
    ///  ### Arguments
    /// * 'day' - Day as an f64
    ///
    /// ### Return
    ///     The coordinates of Earth at the provided time.
    fn ecliptic_cartesian_coords(&self, day: f64) -> Array3d {
        let d = day - 1.5;
        // Julian centuries since J2000.0
        let t = d / 36525.0;
        // Sun's mean longitude, in degrees
        let l_0 = 280.46645 + (36000.76983 * t) + (0.0003032 * t * t);
        // Sun's mean anomaly, in degrees
        let m_0 = 357.52910 + (35999.05030 * t) - (0.0001559 * t * t) - (0.00000048 * t * t * t);

        let c = // Sun's equation of center in degrees
            (1.914600 - 0.004817 * t - 0.000014 * t * t) * sin_deg!(m_0) +
                (0.01993 - 0.000101 * t) * sin_deg!(2f64 * m_0) +
                0.000290 * sin_deg!(3f64 * m_0);

        let ls = l_0 + c; // true elliptical longitude of Sun

        // The eccentricity of the Earth's orbit.
        let e = 0.016708617 - t * (0.000042037 + (0.0000001236 * t));
        // distance from Sun to Earth in astronomical units (AU)
        let distance_in_au = (1.000001018 * (1f64 - e * e)) / (1f64 + e * cos_deg!(m_0 + c));
        let mut x = -distance_in_au * cos_deg!(ls);
        let mut y = -distance_in_au * sin_deg!(ls);

        // Convert AU to Meters
        x *= AU_METER;
        y *= AU_METER;

        // the Earth's center is always on the plane of the ecliptic (z=0), by definition!
        Array3d { x, y, z: 0f64 }
    }

    fn get_solar_object(&self) -> &Solarobj {
        &self.solartype
    }
}

impl KeplerModel for Sun {
    fn ecliptic_cartesian_coords(&self, _day: f64) -> Array3d {
        Array3d {
            x: 0f64,
            y: 0f64,
            z: 0f64,
        }
    }

    fn get_solar_object(&self) -> &Solarobj {
        &self.solartype
    }
}

///  Create the sun.
///
///  ### Return
///       A newly crafted sun object.
fn make_sun() -> Sun {
    let solar_trait = Solarobj::Sun {
        attr: Some(SolarAttr {
            radius: 6.95700e8,
            mass: 1.9891e30,
        }),
    };

    Sun {
        solartype: solar_trait,
        coords: SolarobjCoords::default(),
        velocity: Array3d::default(), // Zero until sim update
    }
}

/// Create the earth.
///
/// ### Argument
/// * 'day' - Day value greater than zero.
///
/// ### Return
///      A newly created earth object.
///
fn make_earth(day: f64) -> Earth {
    let solar_trait = Solarobj::Earth {
        attr: Some(SolarAttr {
            radius: 6.3781e6,
            mass: 5.9722e24,
        }),
    };

    let mut earth_body = Earth {
        solartype: solar_trait,
        coords: SolarobjCoords::default(),
        velocity: Array3d::default(), // Zero until sim update
    };

    let initial_coords = earth_body.ecliptic_cartesian_coords(day);
    earth_body.coords.ahead_coords = initial_coords;
    earth_body.coords.current_coords = initial_coords;
    earth_body.coords.behind_coords = initial_coords;

    earth_body
}

/// Create the moon, geocentric.
///
/// ### Argument
/// * 'day' - Day value greater than zero.
///
/// ### Return
///     A newly created moon PlanetPS object.
///
fn make_moon(day: f64, earth_coords: &Array3d) -> PlanetPS {
    let solar_trait = Solarobj::Moon {
        attr: Some(SolarAttr {
            radius: 1.7381e6,
            mass: 0.07346e24,
        }),
    };

    let mut moon_body = PlanetPS {
        solartype: solar_trait,
        coords: SolarobjCoords::default(),
        velocity: Array3d::default(), // Zero until sim update
        n0: 125.1228,
        nc: -0.0529538083,
        i0: 5.1454,
        ic: 0.0,
        w0: 318.0634,
        wc: 0.1643573223,
        a0: 60.2666 / EARTH_RADII_PER_ASTRONOMICAL_UNIT,
        ac: 0.0,
        e0: 0.054900,
        ec: 0.0,
        m0: 115.3654,
        mc: 13.0649929509,
        mag_base: 0.23,
        mag_phase_factor: 0.026,
        mag_nonlinear_factor: 4.0e-9,
        mag_nonlinear_exponent: 4f64,
    };

    // Calculate the location of the moon and convert to heliocentric coords
    println!("{:?}", day);
    println!("{:?}", moon_body.ecliptic_cartesian_coords(day));
    let initial_coords = moon_body.ecliptic_cartesian_coords(day) + earth_coords;
    moon_body.coords.ahead_coords = initial_coords;
    moon_body.coords.current_coords = initial_coords;
    moon_body.coords.behind_coords = initial_coords;

    moon_body
}
