use crate::{input::SimulationParameters, output};
use crate::types;
use chrono::{DateTime, Duration, TimeZone, Utc};
use serde::{Deserialize, Serialize};
use strum_macros::Display;
use types::Array3d;

const METERS_PER_ASTRONOMICAL_UNIT: f64 = 1.4959787e+11;
const METERS_PER_EARTH_EQUATORIAL_RADIUS: f64 = 6378140.0;
const EARTH_RADII_PER_ASTRONOMICAL_UNIT: f64 =
    METERS_PER_ASTRONOMICAL_UNIT / METERS_PER_EARTH_EQUATORIAL_RADIUS; // 23454.78
const AU_METER: f64 = 1.496e+11;


#[derive(Serialize, Deserialize)]
pub struct InitData {
    pub date: String,                // Datetime in ISO 8601 format
    pub debris: Vec<SimobjT>,         // Debris objects
    pub spacecraft: Vec<SimobjT>, // Spacecraft objects
}

pub trait Simobj {
    fn type_of(&self) -> String;
    fn get_id(&self) -> u32;
    fn id_mut(&mut self) -> &mut u32;
    fn get_ref_coords(&self) -> &types::Array3d;
    fn set_coords(&mut self, value: Array3d);
    fn get_ref_velocity(&self) -> &types::Array3d;
    fn set_velocity(&mut self, value: Array3d);
    fn get_drag_area(&self) -> f64;
    fn get_mass(&self) -> f64;

    fn to_output_form(&self, sim_time: f64) -> output::SimulationObjectParameters {
        let coords = self.get_ref_coords();
        let velocity = self.get_ref_velocity();

        output::SimulationObjectParameters {
            id: self.get_id(),
            sim_time,
            x_coord: coords.x,
            y_coord: coords.y,
            z_coord: coords.z,
            x_velocity: velocity.x,
            y_velocity: velocity.y,
            z_velocity: velocity.z
        }
    }
}

#[derive(Serialize)]
enum SimObjectType{
    Spacecraft,
    Debris
}

impl Default for SimObjectType{
    fn default() -> Self {SimObjectType::Spacecraft}
}

#[derive(Serialize, Deserialize)]
pub struct SimobjT {
    #[serde(skip_deserializing)]
    id: u32,
    #[serde(skip_deserializing)]
    sim_object_type: SimObjectType,
    coords: Array3d,
    velocity: Array3d,
    drag_area: f64,
    mass: f64,
}

impl Simobj for SimobjT {
    fn type_of(&self) -> String {
        match self.sim_object_type{
            Spacecraft => String::from("Spacecraft"),
            Debris => String::from("Debris"),
        }
    }

    fn get_id(&self) -> u32 {
        self.id
    }

    fn id_mut(&mut self) -> &mut u32 {
        &mut self.id
    }


    fn get_ref_coords(&self) -> &Array3d {
        &self.coords
    }

    fn set_coords(&mut self, value: Array3d) {
        self.coords = value;
    }

    fn get_ref_velocity(&self) -> &Array3d {
        &self.velocity
    }

    fn set_velocity(&mut self, value: Array3d) {
        self.velocity = value;
    }

    fn get_drag_area(&self) -> f64 {
        self.drag_area
    }

    fn get_mass(&self) -> f64 {
        self.mass
    }
}


pub struct Environment {
    last_day_update_s: f64,
    sim_time_s: f64, // Simulation time in seconds
    future_day_update_s: f64,
    pub start_time: chrono::DateTime<Utc>,
    // All solar bodies below are synced to future time
    pub sun: Sun,
    last_sun_coords: Array3d, 
    current_sun_coords: Array3d, 
    pub earth: Earth,
    last_earth_coords: Array3d, 
    current_earth_coords: Array3d, 
    pub moon: PlanetPS,
    last_moon_coords: Array3d, 
    current_moon_coords: Array3d, 
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
        *self.sun.mut_coords() = self.sun.ecliptic_cartesian_coords(new_day);
        *self.earth.mut_coords() = self.earth.ecliptic_cartesian_coords(new_day);
        // Calculate new location for moon and convert to heliocentric coords
        *self.moon.mut_coords() = self.moon.ecliptic_cartesian_coords(new_day) + self.earth.get_coords();

    }

    /// Hard update on all solar objects within the simulation. 
    fn update(&mut self, sim_params: &SimulationParameters) {  

        self.last_day_update_s = self.sim_time_s;

        // Calculate postions at current time (postions at future time, TODO optimise)
        let new_time = self.start_time + Duration::seconds(self.sim_time_s as i64);
        self.update_solar_objs(&new_time);


        // Set the corresponding last and current coords
        self.last_sun_coords = self.sun.coords.clone();
        self.current_sun_coords = self.sun.coords.clone();

        self.last_earth_coords = self.earth.coords.clone();
        self.current_earth_coords = self.earth.coords.clone();

        self.last_moon_coords = self.moon.coords.clone();
        self.current_moon_coords = self.moon.coords.clone();


        self.future_day_update_s = self.sim_time_s + (sim_params.sim_solar_step as f64);      

        // Calculate new positions at future time
        new_time = self.start_time + Duration::seconds(self.future_day_update_s as i64);
        self.update_solar_objs(&new_time);

    }

    pub fn sync_sim_objects(&self, sim_objects: &Vec<SimobjT>, sim_params: &SimulationParameters){
        // Perform a linear interpolation from the current solar object coords to the future solar coords
        // Store current solar object coords
        // Set the current solar object coords to the interp value
        // Take stored solar object coords and find delta to current solar object coords
        // Add corrsponding solar object delta to each sim object
        todo!()
    }

    /// Calculates a delta for provided datetime from 0/Jan/2000 00:00 UTC
    ///
    /// ### Argument
    /// * 'datetime' - User provided datetime object.
    ///
    /// ### Return
    ///     The delta from 0/Jan/2000 00:00 UTC in days.
    ///
    fn datetime_to_days(datetime_obj: &DateTime<chrono::Utc>) -> f64 {
        let origin_dt = chrono::Utc.ymd(2000, 1, 1).and_hms(0, 0, 0);

        (1.15741e-5f64 * (*datetime_obj - origin_dt).num_seconds() as f64) as f64
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
        let moon_precalc = make_moon(day, earth_precalc.get_coords());

        let env = Environment {
            last_day_update_s: 0f64,
            sim_time_s: 0f64,
            future_day_update_s: sim_params.sim_solar_step as f64,
            start_time,
            sun: sun_precalc,
            last_sun_coords: sun_precalc.coords,
            current_sun_coords: sun_precalc.coords,
            earth: earth_precalc,
            last_earth_coords: earth_precalc.coords,
            current_earth_coords: earth_precalc.coords,
            moon: moon_precalc,
            last_moon_coords: moon_precalc.coords,
            current_moon_coords: moon_precalc.coords
        };

        env.update(sim_params);

        env
    }
}

#[derive(Display, Clone)]
pub enum Solarobj {
    #[strum(serialize = "sun")]
    Sun { attr: SolarAttr },
    #[strum(serialize = "earth")]
    Earth { attr: SolarAttr },
    #[strum(serialize = "moon")]
    Moon { attr: SolarAttr },
}

#[derive(Clone)]
pub struct SolarAttr {
    radius: f64, // meters
    mass: f64,   // kg
}

impl Solarobj {
    pub fn get_mass_kg(&self) -> f64 {
        match self {
            Solarobj::Sun { attr } => attr.mass,
            Solarobj::Earth { attr } => attr.mass,
            Solarobj::Moon { attr } => attr.mass,
        }
    }
}

pub struct PlanetPS {
    // See  http://www.stjarnhimlen.se/comp/ppcomp.html#4
    solartype: Solarobj, // Type enum of the solar obj
    coords: Array3d,
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

pub struct Earth {
    solartype: Solarobj,
    coords: Array3d,
}

pub struct Sun {
    solartype: Solarobj,
    coords: Array3d,
}

/// Provides utilities for calculating planetary bodies with a Kepler model
mod kepler_utilities {
    use crate::{bodies::{PlanetPS, EARTH_RADII_PER_ASTRONOMICAL_UNIT}, types::Array3d};
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
        Array3d {
            x,
            y,
            z,
        }
    }

    fn distance_between_sim_object(&self, sim_obj: &SimobjT) -> Array3d;

    fn get_coords(&self) -> &Array3d;

    fn mut_coords(&mut self) -> &mut Array3d;

    fn get_solar_object(&self) -> &Solarobj;

    fn to_output_form(&self, sim_time_s: f64) -> output::SolarObjectOut {
        output::SolarObjectOut {
            name: self.get_solar_object().to_string(),
            sim_time: sim_time_s,
            x_coord: self.get_coords().x as f32,
            y_coord: self.get_coords().y as f32,
            z_coord: self.get_coords().z as f32,
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
            _ => Array3d {
                x,
                y,
                z
            },
        }
    }

    fn get_coords(&self) -> &Array3d {
        &self.coords
    }

    fn mut_coords(&mut self) -> &mut Array3d {
        &mut self.coords
    }

    fn get_solar_object(&self) -> &Solarobj {
        &self.solartype
    }

    fn distance_between_sim_object(&self, sim_obj: &SimobjT) -> Array3d {
        &self.coords - sim_obj.get_ref_coords()
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
        Array3d {
            x,
            y,
            z: 0f64,
        }
    }

    fn get_coords(&self) -> &Array3d {
        &self.coords
    }

    fn mut_coords(&mut self) -> &mut Array3d {
        &mut self.coords
    }

    fn get_solar_object(&self) -> &Solarobj {
        &self.solartype
    }

    fn distance_between_sim_object(&self, sim_obj: &SimobjT) -> Array3d {
        &self.coords - sim_obj.get_ref_coords()
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

    fn get_coords(&self) -> &Array3d {
        &self.coords
    }

    fn mut_coords(&mut self) -> &mut Array3d {
        &mut self.coords
    }

    fn get_solar_object(&self) -> &Solarobj {
        &self.solartype
    }

    fn distance_between_sim_object(&self, sim_obj: &SimobjT) -> Array3d {
        &self.coords - sim_obj.get_ref_coords()
    }
}

///  Create the sun.
///
///  ### Return
///       A newly crafted sun object.
fn make_sun() -> Sun {
    let solar_trait = Solarobj::Sun {
        attr: SolarAttr {
            radius: 6.95700e8,
            mass: 1.9891e30,
        },
    };

    Sun {
        solartype: solar_trait,
        coords: Array3d {
            x: 0f64,
            y: 0f64,
            z: 0f64,
        },
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
        attr: SolarAttr {
            radius: 6.3781e6,
            mass: 5.9722e24,
        },
    };

    let mut earth_body = Earth {
        solartype: solar_trait,
        coords: Array3d {
            x: 0f64,
            y: 0f64,
            z: 0f64,
        },
    };

    earth_body.coords = earth_body.ecliptic_cartesian_coords(day);

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
        attr: SolarAttr {
            radius: 1738.1,
            mass: 0.07346e24,
        },
    };

    let mut moon_body = PlanetPS {
        solartype: solar_trait,
        coords: Array3d {
            x: 0f64,
            y: 0f64,
            z: 0f64,
        },
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
    moon_body.coords = moon_body.ecliptic_cartesian_coords(day) + earth_coords;

    moon_body
}
