use std::f64::consts::PI;

use chrono::{DateTime, Datelike, Timelike, Utc};
use serde::{Deserialize, Serialize};
use strum_macros::Display;
use types::Array3d;

use crate::{input::SwIndex, types::LLH};

use super::common::ct2lst;

const METERS_PER_ASTRONOMICAL_UNIT: f64 = 1.4959787e+11;
const METERS_PER_EARTH_EQUATORIAL_RADIUS: f64 = 6378140.0;
const EARTH_RADII_PER_ASTRONOMICAL_UNIT: f64 =
    METERS_PER_ASTRONOMICAL_UNIT / METERS_PER_EARTH_EQUATORIAL_RADIUS; // 23454.78
const AU_METER: f64 = 1.496e+11;

#[derive(Display, Clone, Deserialize, Serialize)]
#[serde(tag = "type")]
pub enum Solarobj {
    #[strum(serialize = "sun")]
    Sun,
    #[strum(serialize = "earth")]
    Earth,
    #[strum(serialize = "moon")]
    Moon,
}

#[derive(Clone, Deserialize, Serialize)]
pub struct SolarAttr {
    pub eqradius: f64,  // meters
    pub mass: f64,      // kg
    pub obliquity: f64, // degrees relative to the ecliptic
}

#[derive(Clone)]
pub struct ModelState {
    pub coords: SolarobjCoords,
    pub velocity: Array3d,
}

#[derive(Clone)]
#[allow(dead_code)]
pub struct PlanetPSModel {
    pub state: ModelState,
    solartype: Solarobj,
    // See  http://www.stjarnhimlen.se/comp/ppcomp.html#4
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
pub struct EarthModel {
    pub state: ModelState,
}

#[derive(Clone)]
pub struct SunModel {
    pub state: ModelState,
}

pub struct Sun {
    pub model: SunModel,
    pub attr: &'static SolarAttr,
}

pub struct Earth {
    sw_indices: Vec<SwIndex>, // Space Weather Indices.
    current_sw: usize,
    pub model: EarthModel,
    pub attr: &'static SolarAttr,
}

pub struct Moon {
    pub model: PlanetPSModel,
    pub attr: &'static SolarAttr,
}

/// Provides utilities for calculating planetary bodies with a Kepler model.
/// Most of the code used to calculate planetary body coordinates was referenced
/// from source for http://cosinekitty.com/solar_system.html. See, http://cosinekitty.com/astronomy.js.
mod kepler_utilities {
    use crate::types::Array3d;
    use std::f64::{self, consts};

    use super::{PlanetPSModel, EARTH_RADII_PER_ASTRONOMICAL_UNIT};

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
        356.0470 + (0.9856002585 * day)
    }

    /// Calculates the argument of perihelion for the Sun.
    fn sun_argument_of_perihelion(day: f64) -> f64 {
        282.9404 + (4.70935e-5 * day)
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

    pub fn lunar_pertub(body: &PlanetPSModel, xh: f64, yh: f64, zh: f64, day: f64) -> Array3d {
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
            x: xp,
            y: yp,
            z: zp,
        }
    }
}

pub trait KeplerModel {
    fn ecliptic_cartesian_coords(&self, day: f64) -> Array3d;

    fn perturb(&self, x: f64, y: f64, z: f64, _day: f64) -> Array3d {
        Array3d { x, y, z }
    }
}

#[derive(Clone, Copy, Default)]
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
    pub fn lerp_set(&mut self, interp_point: f64) {
        // Precise method, which guarantees v = v1 when t = 1.
        self.current_coords =
            ((1f64 - interp_point) * self.behind_coords) + interp_point * self.ahead_coords;
    }

    /// Calculate the velocity of the solar obj from behind to ahead coords.
    ///
    /// ### Arguments
    /// * 'step_time_s' - Sim solar step time in seconds
    ///
    /// ### Returns
    /// The velocity of the solar obj
    ///
    pub fn calc_velocity_full_range(&self, step_time_s: f64) -> Array3d {
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
    pub fn calc_velocity_relative_full_range(
        &self,
        relative: &SolarobjCoords,
        step_time_s: f64,
    ) -> Array3d {
        ((self.ahead_coords - relative.ahead_coords)
            - (self.behind_coords - relative.behind_coords))
            / step_time_s
    }
}

impl KeplerModel for PlanetPSModel {
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

        self.perturb(xh, yh, zh, day)
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
            Solarobj::Moon => kepler_utilities::lunar_pertub(self, x, y, z, day),
            _ => Array3d { x, y, z },
        }
    }
}

impl PlanetPSModel {
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

impl KeplerModel for EarthModel {
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
}

impl KeplerModel for SunModel {
    fn ecliptic_cartesian_coords(&self, _day: f64) -> Array3d {
        Array3d {
            x: 0f64,
            y: 0f64,
            z: 0f64,
        }
    }
}

impl Sun {
    pub fn new() -> Self {
        Sun {
            model: SunModel {
                state: ModelState {
                    coords: SolarobjCoords::default(),
                    velocity: Array3d::default(),
                },
            },
            attr: &SolarAttr {
                eqradius: 6.95700e8,
                mass: 1.9891e30,
                obliquity: 0.0,
            },
        }
    }
}

impl Earth {
    /// Grab Earth space weather data for use in earth based models.
    ///
    ///  ### Arguments
    /// * 'query_time' - UTC Datetime object used to find space weather index.
    ///
    /// ### Return
    ///     Space weather index matching provided query time if found.
    ///
    /// ### Panics
    ///     Panics if a match cannot be found.
    ///
    fn get_space_weather_index(&self, query_time: DateTime<chrono::Utc>) -> usize {
        let cmp_func = |probe: &SwIndex| probe.cmp(query_time);
        let search_result = self.sw_indices.binary_search_by(cmp_func);

        // A match should always be found if correct simulation configuration has been provided.
        match search_result {
            Ok(index) => index,
            Err(_) => panic!(),
        }
    }

    pub fn set_query_space_weather_index(&mut self, query_time: DateTime<chrono::Utc>) {
        self.current_sw = self.get_space_weather_index(query_time);
    }

    /// Call into the nrlmsise00 model for generating atmospheric parameters.
    /// This method relies on Earth state being set correctly before calling.
    ///
    ///  ### Arguments
    /// * 'current_datetime' - UTC Datetime of to use within the model call.
    /// * 'fixed_coords' - ECEF coordinates of query point.
    ///
    /// ### Return
    ///     Output object containing data generated from the nrlmsise00 model.
    ///
    pub fn nrlmsise00_model(
        &self,
        current_datetime: DateTime<Utc>,
        fixed_coords: &LLH,
    ) -> nrlmsise00c::NRLMSISEOutput {
        let seconds_from_midnight = current_datetime.num_seconds_from_midnight() as f64;
        let sw_index = self.sw_indices.get(self.current_sw).unwrap();

        let mut input = nrlmsise00c::NRLMSISEInput {
            year: current_datetime.year(),
            doy: current_datetime.ordinal() as i32,
            sec: seconds_from_midnight,
            alt: fixed_coords.alt / 1000.0,
            g_lat: fixed_coords.lat,
            g_long: fixed_coords.long,
            lst: seconds_from_midnight / 3600.0 + fixed_coords.long / 15.0, // (lst=sec/3600 + g_long/15)
            f107A: sw_index.4,
            f107: sw_index.3,
            ap: sw_index.2,
            ap_a: [0f64; 7usize],
        };
        let flags = nrlmsise00c::NRLMSISEFlags {
            switches: [1; 24usize],
            sw: [0f64; 24usize],
            swc: [0f64; 24usize],
        };

        nrlmsise00c::gtd7_safe(&mut input, &flags)
    }

    /// Convert Earth-centered inertial coordinates to geographic spherical coords
    /// Ref: https://idlastro.gsfc.nasa.gov/ftp/pro/astro/eci2geo.pro
    ///
    /// ### Arguments
    /// * 'eci_coord' - ECI coordinate to be translated into Earth fixed coordinates.
    ///                 Coordinates should be provided in meters.
    /// * 'jd2000' - Julian day in decimal hours from the J2000 epoch.
    ///
    /// ### Returns
    ///     LLH struct containing latitude, longitude and altitude (in meters).
    ///  
    pub fn eci2geo(&self, eci_coord: &Array3d, jd2000: f64) -> LLH {
        let eci_coord_km = eci_coord / 1000.0;

        let re = self.attr.eqradius / 1000.0;
        let theta = eci_coord_km.y.atan2(eci_coord_km.x); // azimuth
        let gst = ct2lst(0.0, jd2000); // Greenwich mean sidereal time

        let angle_sid = gst * 2.0 * PI / 24.0; // sidereal angle
        let long = (theta - angle_sid).rem_euclid(2.0 * PI).to_degrees(); // longitude
        let r = f64::sqrt(eci_coord_km.x.powf(2.0) + eci_coord_km.y.powf(2.0));
        let lat = eci_coord_km.z.atan2(r); // latitude
        let alt = r / lat.cos() - re; // altitude

        LLH {
            lat: lat.to_degrees(),
            long,
            alt: alt * 1000.0, // Convert to meters to keep consistency.
        }
    }

    pub fn new(day: f64, sw_indices: &Vec<SwIndex>) -> Self {
        let mut sw_indices_clone = sw_indices.clone();
        // Ensure the indices list is sorted as binary search will be done on this data.
        sw_indices_clone.sort_by_key(|k| k.0);

        let mut earth_body = Self {
            model: EarthModel {
                state: ModelState {
                    coords: SolarobjCoords::default(),
                    velocity: Array3d::default(),
                },
            },
            sw_indices: sw_indices_clone,
            current_sw: 0,
            attr: &SolarAttr {
                eqradius: 6378137.0,
                mass: 5.9722e24,
                obliquity: 23.44,
            },
        };

        let initial_coords = earth_body.model.ecliptic_cartesian_coords(day);
        earth_body.model.state.coords.ahead_coords = initial_coords;
        earth_body.model.state.coords.current_coords = initial_coords;
        earth_body.model.state.coords.behind_coords = initial_coords;

        earth_body
    }
}

impl Moon {
    pub fn new(day: f64, earth_coords: &Array3d) -> Self {
        let mut moon_body = Self {
            model: PlanetPSModel {
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
                state: ModelState {
                    coords: SolarobjCoords::default(),
                    velocity: Array3d::default(), // Zero until sim update
                },
                solartype: Solarobj::Moon,
            },
            attr: &SolarAttr {
                eqradius: 1.7381e6,
                mass: 7.3459e22,
                obliquity: 1.5424,
            },
        };

        // Calculate the location of the moon and convert to heliocentric coords
        let initial_coords = moon_body.model.ecliptic_cartesian_coords(day) + earth_coords;
        moon_body.model.state.coords.ahead_coords = initial_coords;
        moon_body.model.state.coords.current_coords = initial_coords;
        moon_body.model.state.coords.behind_coords = initial_coords;

        moon_body
    }
}
