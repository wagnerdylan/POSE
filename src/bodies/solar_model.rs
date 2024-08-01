use std::f64::consts::PI;

use chrono::{DateTime, Datelike, Timelike, Utc};
use serde::{Deserialize, Serialize};
use strum_macros::Display;
use types::Array3d;

use crate::{input::SwIndex, types::LLH};

use super::common::ct2lst;

pub fn furnish_spice() {
    spice::furnsh("data/spice/spk/planets/de440s.bsp");
    spice::furnsh("data/spice/lsk/latest_leapseconds.tls");
    spice::furnsh("data/spice/pck/earth_200101_990825_predict.bpc");
}

#[derive(Display, Clone, Deserialize, Serialize, Default, PartialEq)]
#[serde(tag = "type")]
pub enum Solarobj {
    #[strum(serialize = "sun")]
    Sun,
    #[default]
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
    pub coords: Array3d,
    pub velocity: Array3d,
}

#[derive(Clone)]
pub struct SpiceSunModel {
    pub state: ModelState,
}

#[derive(Clone)]
pub struct Sun {
    pub model: SpiceSunModel,
    pub attr: &'static SolarAttr,
}

#[derive(Clone)]
pub struct SpiceEarthModel {
    pub state: ModelState,
}

#[derive(Clone)]
pub struct Earth {
    sw_indices: Vec<SwIndex>, // Space Weather Indices.
    current_sw: usize,
    pub model: SpiceEarthModel,
    pub attr: &'static SolarAttr,
}

#[derive(Clone)]
pub struct SpiceMoonModel {
    pub state: ModelState,
}

#[derive(Clone)]
pub struct Moon {
    pub model: SpiceMoonModel,
    pub attr: &'static SolarAttr,
}

pub trait OrbitModel {
    /// Compute positions X, Y and Z along with velocity values VX, VY and VZ.
    ///
    /// ### Arguments
    /// * 'et' - Epoch time in seconds from J2000.
    ///
    /// ### Return
    ///     A tuple of Array3d objects containing position and velocity respectfully.
    ///
    fn ecliptic_cartesian_coords(&self, et: f64) -> (Array3d, Array3d);

    /// Update internal model state within the object following the epoch time
    /// parameter.
    ///
    /// ### Arguments
    /// * 'et' - Epoch time in seconds from J2000.
    ///
    fn update_position_velocity(&mut self, et: f64);
}

/// Convert three values of a given slice to an Array3d object and
/// scale the output from km to m to match what POSE uses as unit scaling.
///
/// ### Arguments
/// * 'value' - Slice of three values to be packed and scaled.
///
/// ### Return
///     The final scaled and packed Array3d object.
///
fn scale_and_pack_spk(value: &[f64]) -> Array3d {
    const KM_TO_M: f64 = 1000.0;
    Array3d {
        x: *value.get(0).unwrap() * KM_TO_M,
        y: *value.get(1).unwrap() * KM_TO_M,
        z: *value.get(2).unwrap() * KM_TO_M,
    }
}

/// Calculate position and velocity for a specified solar body in the
/// J2000 frame using the Sun as the observer.
///
/// ### Arguments
/// * 'target' - String containing the name of the body in which to calculate position and velocity.
///              This string must be of a variant in which SPICE may accept.
/// * 'et' - Epoch time in seconds from J2000.
///
/// ### Return
///     A tuple of Array3d objects containing position and velocity respectfully.
///
fn spk_j2000_sun_ezr(target: &str, et: f64) -> (Array3d, Array3d) {
    let (output, _) = spice::spkezr(target, et, "J2000", "NONE", "SUN");
    (
        // The first three elements of the spkezr function output are cartesian coordinates X, Y and Z.
        scale_and_pack_spk(&output[0..3]),
        // The last three elements of the spkezr function output are velocity values VX, VY and VZ.
        scale_and_pack_spk(&output[3..6]),
    )
}

impl OrbitModel for SpiceSunModel {
    fn ecliptic_cartesian_coords(&self, et: f64) -> (Array3d, Array3d) {
        spk_j2000_sun_ezr("SUN", et)
    }

    fn update_position_velocity(&mut self, et: f64) {
        let (pos, vel) = self.ecliptic_cartesian_coords(et);
        self.state.coords = pos;
        self.state.velocity = vel;
    }
}

impl OrbitModel for SpiceEarthModel {
    fn ecliptic_cartesian_coords(&self, et: f64) -> (Array3d, Array3d) {
        spk_j2000_sun_ezr("EARTH", et)
    }

    fn update_position_velocity(&mut self, et: f64) {
        let (pos, vel) = self.ecliptic_cartesian_coords(et);
        self.state.coords = pos;
        self.state.velocity = vel;
    }
}

impl OrbitModel for SpiceMoonModel {
    fn ecliptic_cartesian_coords(&self, et: f64) -> (Array3d, Array3d) {
        spk_j2000_sun_ezr("MOON", et)
    }

    fn update_position_velocity(&mut self, et: f64) {
        let (pos, vel) = self.ecliptic_cartesian_coords(et);
        self.state.coords = pos;
        self.state.velocity = vel;
    }
}

impl Sun {
    pub fn new(et: f64) -> Self {
        let mut sun_body = Sun {
            model: SpiceSunModel {
                state: ModelState {
                    coords: Array3d::default(),
                    velocity: Array3d::default(),
                },
            },
            attr: &SolarAttr {
                eqradius: 6.95700e8,
                mass: 1.9891e30,
                obliquity: 0.0,
            },
        };
        sun_body.model.update_position_velocity(et);

        sun_body
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

    pub fn new(et: f64, sw_indices: &Vec<SwIndex>) -> Self {
        let mut sw_indices_clone = sw_indices.clone();
        // Ensure the indices list is sorted as binary search will be done on this data.
        sw_indices_clone.sort_by_key(|k| k.0);

        let mut earth_body = Self {
            model: SpiceEarthModel {
                state: ModelState {
                    coords: Array3d::default(),
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
        earth_body.model.update_position_velocity(et);

        earth_body
    }
}

impl Moon {
    pub fn new(et: f64) -> Self {
        let mut moon_body = Self {
            model: SpiceMoonModel {
                state: ModelState {
                    coords: Array3d::default(),
                    velocity: Array3d::default(),
                },
            },
            attr: &SolarAttr {
                eqradius: 1.7381e6,
                mass: 7.3459e22,
                obliquity: 1.5424,
            },
        };
        moon_body.model.update_position_velocity(et);

        moon_body
    }
}
