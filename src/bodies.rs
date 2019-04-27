use serde::{Deserialize, Serialize};


const METERS_PER_ASTRONOMICAL_UNIT: f32 = 1.4959787e+11;
const METERS_PER_EARTH_EQUATORIAL_RADIUS: f32 = 6378140.0;
const EARTH_RADII_PER_ASTRONOMICAL_UNIT: f32 =
    METERS_PER_ASTRONOMICAL_UNIT / METERS_PER_EARTH_EQUATORIAL_RADIUS;      // 23454.78


pub type SimobjT = Box<dyn Simobj>;
pub type PlanetBody = Box<dyn KeplerModel>;

#[derive(Serialize, Deserialize, Clone)]
pub struct InitData {
    pub date: String,
    pub debris: Vec<Debris>,
    pub spacecraft: Vec<Spacecraft>
}

pub trait Simobj {
    fn type_of(&self) -> String;
    fn get_id(&self) -> u32;
    fn id_mut(&mut self) -> &mut u32;
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Spacecraft{
    #[serde(skip_deserializing)]
    id: u32,
    x_dis: f64,
    y_dis: f64,
    z_dis: f64,
    x_vel: f64,
    y_vel: f64,
    z_vel: f64
    // Add more
}

impl Simobj for Spacecraft {
    fn type_of(&self) -> String {
        return String::from("Spacecraft");
    }
    fn get_id(&self) -> u32 { self.id }

    fn id_mut(&mut self) -> &mut u32 {&mut self.id}
}

/// Struct for holding attributes relating to debris
#[derive(Serialize, Deserialize, Clone)]
pub struct Debris{
    #[serde(skip_deserializing)]
    id: u32,
    x_dis: f64,
    y_dis: f64,
    z_dis: f64,
    x_vel: f64,
    y_vel: f64,
    z_vel: f64
    // Add more
}

impl Simobj for Debris {
    fn type_of(&self) -> String {
        return String::from("Debris");
    }

    fn get_id(&self) -> u32 { self.id }

    fn id_mut(&mut self) -> &mut u32 {&mut self.id}
}

#[derive(Debug)]
pub enum Solarobj{
    Sun{attr: SolarAttr},
    Earth{attr: SolarAttr},
    Moon{attr: SolarAttr}
}

#[derive(Debug)]
pub struct SolarAttr{
    radius: f64, // meters
    mass: f64 // kg
}

#[derive(Debug)]
pub struct PlanetPS { // See  http://www.stjarnhimlen.se/comp/ppcomp.html#4
    solartype: Solarobj, // Type enum of the solar obj
    coords: CartesianCoords,
    n0: f32, nc: f32, // N0 = longitude of the ascending node (deg).  Nc = rate of change in deg/day
    i0: f32, ic: f32, // inclination to the ecliptic (deg)
    w0: f32, wc: f32, // argument of perihelion (deg)
    a0: f32, ac: f32, // semi-major axis, or mean distance from Sun (AU)
    e0: f32, ec: f32, // eccentricity (0=circle, 0..1=ellipse, 1=parabola)
    m0: f32, mc: f32, // M0 = mean anomaly  (deg) (0 at perihelion; increases uniformly with time).  Mc ("mean motion") = rate of change
    mag_base: f32,
    mag_phase_factor: f32,
    mag_nonlinear_factor: f32,
    mag_nonlinear_exponent: f32
}

#[derive(Debug)]
pub struct Earth {
    solartype: Solarobj,
    coords: CartesianCoords
}

#[derive(Debug)]
pub struct Sun {
    solartype: Solarobj,
    coords: CartesianCoords
}

#[derive(Debug)]
pub struct CartesianCoords {
    is_meters: bool,
    heliocentric: bool, // False if geocentric
    pub xh: f32, // X location in meters
    pub yh: f32, // Y location in meters
    pub zh: f32  // Z location in meters
}

impl CartesianCoords {
    
    /// Converts the cartesian coords from Au to meters.
    pub fn to_meters(&mut self){
        if !self.is_meters {
            self.is_meters = true;
            self.xh *= 149600000000f32;
            self.zh *= 149600000000f32;
            self.yh *= 149600000000f32;
        }
    }

    ///Converts the cartesian coords from meters to AU.
    pub fn to_au(&mut self){
        if self.is_meters {
            self.is_meters = false;
            self.xh /= 149600000000f32;
            self.zh /= 149600000000f32;
            self.yh /= 149600000000f32;
        }
    }
}

/// Provides utilities for calculating planetary bodies with a Kepler model
mod kepler_utilities {
    use std::f32::{self, consts};
    use crate::bodies::{PlanetPS, KeplerModel, CartesianCoords, EARTH_RADII_PER_ASTRONOMICAL_UNIT};


    /// Calculate the eccentric anomaly for a given body.
    /// ### Arguments
    /// * 'e' - TODO
    /// * 'm' - TODO
    /// 
    /// ### Returns
    ///      The eccentric anomaly for the provided input parameters.
    pub fn eccentric_anomaly(e: f32, m: f32) -> f32 {

        let deg_from_rad = 180f32 / consts::PI;
        let mut ecc: f32 = m + (e * sin_deg!(m) * (1f32 + (e * cos_deg!(m))));

        loop {
            let f: f32 = ecc - (ecc - (deg_from_rad * e * sin_deg!(ecc)) - m) / (1f32 - e * cos_deg!(ecc));
            let error = (f - ecc).abs();
            ecc = f;

            if error < 1.0e-2 { break; }
        };

        ecc
    }
    
    /// Calculates the mean anomaly for the Sun.
    fn mean_anomaly_of_sun(day: f32) -> f32 {
        356.0470 + (0.9856002585 * day)
    }


    /// Calculates the argument of perihelion for the Sun.
    fn sun_argument_of_perihelion(day: f32) -> f32 {
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
    fn ecliptic_lat_lon(xh: f32, yh: f32, zh: f32) -> (f32, f32){
        (atan2_deg!(yh, xh), atan2_deg!(zh, (xh*xh + yh*yh).sqrt()))
    }

    pub fn lunar_pertub(body: &PlanetPS, xh: f32, yh: f32, zh: f32, day: f32) -> CartesianCoords{
        let ms = mean_anomaly_of_sun(day); // mean anomaly of Sun
        let ws = sun_argument_of_perihelion(day); // Sun's argument of perihelion
        let ls = ms + ws; // mean longitude of Sun

        let mm = body.mean_anomaly(day); // Moon's mean anomaly
        let nm = body.node_longitude(day); // longitude of Moon's node
        let wm = body.perihelion(day); // Moon's argument of perihelion
        let lm = mm + wm + nm; // Mean longitude of the Moon

        let d = lm - ls; // mean elongation of the Moon
        let f = lm - nm; // argument of latitude for the Moon

        let delta_long =
           -1.274 * sin_deg!(mm - 2f32*d)       +   // the Evection
            0.658 * sin_deg!(2f32*d)            -   // the Variation
            0.186 * sin_deg!(ms)                -   // the Yearly Equation
            0.059 * sin_deg!(2f32*mm - 2f32*d)  -
            0.057 * sin_deg!(mm - 2f32*d + ms)  +
            0.053 * sin_deg!(mm + 2f32*d)       +
            0.046 * sin_deg!(2f32*d - ms)       +
            0.041 * sin_deg!(mm - ms)           -
            0.035 * sin_deg!(d)                 -   // the Parallactic Equation
            0.031 * sin_deg!(mm + ms)           -
            0.015 * sin_deg!(2f32*f - 2f32*d)   +
            0.011 * sin_deg!(mm - 4f32*d)
        ;

        let delta_lat =
           -0.173 * sin_deg!(f - 2f32*d)        -
            0.055 * sin_deg!(mm - f - 2f32*d)   -
            0.046 * sin_deg!(mm + f - 2f32*d)   +
            0.033 * sin_deg!(f + 2f32*d)        +
            0.017 * sin_deg!(2f32*mm + f)
        ;

        let delta_radius =
           -0.58 * cos_deg! (mm - 2f32*d)       -
            0.46 * cos_deg! (2f32*d)
        ;

        let (mut lonecl, mut latecl) = ecliptic_lat_lon(xh, yh, zh);

        let mut r = (xh*xh + yh*yh + zh*zh).sqrt();

        lonecl += delta_long;
        latecl += delta_lat;
        r += delta_radius / EARTH_RADII_PER_ASTRONOMICAL_UNIT;

        let coslon = cos_deg! (lonecl);
        let sinlon = sin_deg! (lonecl);
        let coslat = cos_deg! (latecl);
        let sinlat = sin_deg! (latecl);

        let xp = r * coslon * coslat;
        let yp = r * sinlon * coslat;
        let zp = r * sinlat;

        CartesianCoords{xh: xp, yh: yp, zh: zp, is_meters: false, heliocentric: false}
    }

}

pub trait KeplerModel{

    fn ecliptic_cartesian_coords(&self, day: f32) -> CartesianCoords;

    fn perturb(&self, xh: f32, yh: f32, zh: f32, day: f32) -> CartesianCoords {
        CartesianCoords{xh, yh, zh, is_meters: false, heliocentric: true}
    }

    fn mean_anomaly(&self, day: f32) -> f32;

    fn node_longitude(&self, day: f32) -> f32;

    fn perihelion(&self, day: f32) -> f32;

    fn get_coords(&self) -> &CartesianCoords;

    fn mut_coords(&mut self) -> &mut CartesianCoords;

}

impl KeplerModel for PlanetPS{

    fn ecliptic_cartesian_coords(&self, day: f32) -> CartesianCoords {
        // Default impl
        let a = self.a0 + (day * self.ac);
        let e = self.e0 + (day * self.ec);
        let m_u = self.m0 + (day * self.mc);
        let n_u = self.n0 + (day * self.nc);
        let w = self.w0 + (day * self.wc);
        let i = self.i0 + (day * self.ic);
        let ecc = kepler_utilities::eccentric_anomaly(e, m_u);

        let xv = a * (cos_deg!(ecc) - e);
        let yv = a * ((1.0f32 - e*e).sqrt() * sin_deg!(ecc));

        let v = atan2_deg!(yv, xv); // True anomaly in degrees: the angle from perihelion of the body as seen by the Sun.
        let r = (xv*xv + yv*yv).sqrt(); // Distance from the Sun to the planet in AU

        let cos_n  = cos_deg!(n_u);
        let sin_n  = sin_deg!(n_u);
        let cosi  = cos_deg!(i);
        let sini  = sin_deg!(i);
        let cos_vw = cos_deg!(v + w);
        let sin_vw = sin_deg!(v + w);
    
        // Now we are ready to calculate (unperturbed) ecliptic cartesian heliocentric coordinates.
        let xh = r * (cos_n * cos_vw - sin_n * sin_vw * cosi);
        let yh = r * (sin_n * cos_vw + cos_n * sin_vw * cosi);
        let zh = r * sin_vw * sini;

        let mut coords = self.perturb(xh, yh, zh, day);
        coords
    }
    
    /// Calculates additional perturbations on top of main heliocentric position calculation.
    /// Matches PlanetPS bodies using the type enum.
    ///  
    /// ### Arguments
    ///  * 'xh' - X coord
    ///  * 'yh' - Y coord
    ///  * 'zh' - Z coord
    ///  8 'day' - Day value
    /// 
    /// ### Returns
    ///      Cartesian coords with the added perturbations.
    /// 
    fn perturb(&self, xh: f32, yh: f32, zh: f32, day: f32) -> CartesianCoords {
        match &self.solartype {
            Solarobj::Moon{attr} => kepler_utilities::lunar_pertub(self, xh, yh, zh, day),
            _ => CartesianCoords{xh, yh, zh, is_meters: false, heliocentric: true}
        }
    }

    fn mean_anomaly(&self, day: f32) -> f32 {
        self.m0 + (day * self.mc)
    }

    fn node_longitude(&self, day: f32) -> f32 {
        self.n0 + (day * self.nc)
    }

    fn perihelion(&self, day: f32) -> f32 {
        self.w0 + (day * self.wc)
    }

    fn get_coords(&self) -> &CartesianCoords {
        &self.coords
    }

    fn mut_coords(&mut self) -> &mut CartesianCoords {
        &mut self.coords
    }

}

impl KeplerModel for Earth {
    
    /// Calculate the position of Earth relative to the Sun.
    /// Calls function earth_ecliptic_cartesian_coords in kepler_utilities
    /// 
    ///  ### Arguments
    /// * 'day' - Day as an f32
    /// 
    /// ### Return
    ///     The coordinates of Earth at the provided time.
    fn ecliptic_cartesian_coords(&self, day: f32) -> CartesianCoords {

        let d = day - 1.5;
        // Julian centuries since J2000.0
        let t = d / 36525.0;
        // Sun's mean longitude, in degrees
        let l_0 = 280.46645 + (36000.76983 * t) + (0.0003032 * t * t);
        // Sun's mean anomaly, in degrees
        let m_0 = 357.52910 + (35999.05030 * t) - (0.0001559 * t * t) - (0.00000048 * t * t * t);

        let c = // Sun's equation of center in degrees
            (1.914600 - 0.004817 * t - 0.000014 * t * t) * sin_deg!(m_0) +
                (0.01993 - 0.000101 * t) * sin_deg!(2f32 * m_0) +
                0.000290 * sin_deg!(3f32 * m_0);

        let ls = l_0 + c; // true elliptical longitude of Sun

        // The eccentricity of the Earth's orbit.
        let e = 0.016708617 - t * (0.000042037 + (0.0000001236 * t));
        // distance from Sun to Earth in astronomical units (AU)
        let distance_in_au = (1.000001018 * (1f32 - e * e)) / (1f32 + e * cos_deg!(m_0 + c));
        let x = -distance_in_au * cos_deg!(ls);
        let y = -distance_in_au * sin_deg!(ls);

        // the Earth's center is always on the plane of the ecliptic (z=0), by definition!
        let coords = CartesianCoords { xh: x, yh: y, zh: 0f32, is_meters: false, heliocentric: true };

        coords
    }

    fn mean_anomaly(&self, day: f32) -> f32 {
        0f32
    }

    fn node_longitude(&self, day: f32) -> f32 {
        0f32
    }

    fn perihelion(&self, day: f32) -> f32 {
        0f32
    }

    fn get_coords(&self) -> &CartesianCoords {
        &self.coords
    }

    fn mut_coords(&mut self) -> &mut CartesianCoords {
        &mut self.coords
    }
}


impl KeplerModel for Sun {

    fn ecliptic_cartesian_coords(&self, day:f32) -> CartesianCoords {
        CartesianCoords{xh: 0f32, yh: 0f32, zh: 0f32, is_meters: false, heliocentric: true}
    }

    fn mean_anomaly(&self, day: f32) -> f32 {
        0f32
    }

    fn node_longitude(&self, day: f32) -> f32 {
        0f32
    }

    fn perihelion(&self, day: f32) -> f32 {
        0f32
    }

    fn get_coords(&self) -> &CartesianCoords {
        &self.coords
    }

    fn mut_coords(&mut self) -> &mut CartesianCoords {
        &mut self.coords
    }
}

///  Create the sun.
/// 
///  ### Return
///       A newly crafted sun object.
fn make_sun() -> Sun {

    let solar_trait = Solarobj::Sun{attr: SolarAttr{radius: 6.95700e8, mass: 1.9891e30}};

    let sun_body = Sun{solartype: solar_trait, coords: CartesianCoords{xh: 0f32, yh: 0f32, zh: 0f32, is_meters: false, heliocentric: true}};

    sun_body
}

/// Create the earth.
///
/// ### Argument
/// * 'day' - Day value greater than zero.
///
/// ### Return
///      A newly created earth object.
///
fn make_earth(day: f32) -> Earth {

    let solar_trait = Solarobj::Earth{attr: SolarAttr{radius: 6.3781e6, mass: 5.9722e24}};

    let mut earth_body = Earth{solartype: solar_trait,
                           coords: CartesianCoords{xh: 0f32, yh: 0f32, zh: 0f32, 
                                                   is_meters: false, heliocentric: true}};

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
fn make_moon(day: f32) -> PlanetPS {

    let solar_trait = Solarobj::Moon {attr: SolarAttr{radius: 1738.1, mass: 0.07346e24}};

    let mut moon_body = PlanetPS{
        solartype: solar_trait,
        coords: CartesianCoords{xh: 0f32, yh: 0f32, zh: 0f32, is_meters: false, heliocentric: false},
        n0: 125.1228,   nc: -0.0529538083,
        i0: 5.1454,     ic: 0.0,
        w0: 318.0634,   wc: 0.1643573223,
        a0: 60.2666 / EARTH_RADII_PER_ASTRONOMICAL_UNIT,
        ac: 0.0,
        e0: 0.054900,   ec: 0.0,
        m0: 115.3654,   mc: 13.0649929509,
        mag_base: 0.23,
        mag_phase_factor: 0.026,
        mag_nonlinear_factor: 4.0e-9,
        mag_nonlinear_exponent: 4f32
    };

    moon_body.coords = moon_body.ecliptic_cartesian_coords(day);

    moon_body
}

/// Creates the initial vector of solar system objects.
/// 0 - Sun, 1 - Earth, 2 - Moon
///
/// ### Argument
/// * 'day' - The day value greater than zero. From 2000-01-01
///
/// ### Return
///      A vector with the elements defined above.
///
pub fn solar_system_objs(day: f32) -> Vec<PlanetBody> {
    let mut solar_bodies: Vec<PlanetBody> = Vec::new();

    solar_bodies.push(Box::new(make_sun()));
    solar_bodies.push(Box::new(make_earth(day)));
    solar_bodies.push(Box::new(make_moon(day)));

    solar_bodies
}

/// Updates the coords for all PlanetBody objects in the provided vector.
///
/// ### Argument
/// * 'ss_objs' - PlanetBody objects to be updated.
///
pub fn update_solar_system_objs(ss_objs: &mut Vec<PlanetBody>, day: f32){

    for obj in ss_objs {
        *obj.mut_coords() = obj.ecliptic_cartesian_coords(day);
    }
}
