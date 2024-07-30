use chrono::{DateTime, Datelike, Timelike, Utc};

use crate::types::Array3d;

pub fn equatorial_to_ecliptic(x: &Array3d, obliquity: f64) -> Array3d {
    let r1 = Array3d {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    let r2 = Array3d {
        x: 0.0,
        y: cos_deg!(obliquity),
        z: sin_deg!(obliquity),
    };
    let r3 = Array3d {
        x: 0.0,
        y: -sin_deg!(obliquity),
        z: cos_deg!(obliquity),
    };

    Array3d {
        x: x.dot(&r1),
        y: x.dot(&r2),
        z: x.dot(&r3),
    }
}

pub fn ecliptic_to_equatorial(x: &Array3d, obliquity: f64) -> Array3d {
    let r1 = Array3d {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    let r2 = Array3d {
        x: 0.0,
        y: cos_deg!(obliquity),
        z: -sin_deg!(obliquity),
    };
    let r3 = Array3d {
        x: 0.0,
        y: sin_deg!(obliquity),
        z: cos_deg!(obliquity),
    };

    Array3d {
        x: x.dot(&r1),
        y: x.dot(&r2),
        z: x.dot(&r3),
    }
}

/// Convert UTC datetime objects into julian days.
/// Ref: https://idlastro.gsfc.nasa.gov/ftp/pro/astro/jdcnv.pro
#[inline]
fn jdconv(date: &DateTime<Utc>) -> f64 {
    let l: i32 = (date.month() as i32 - 14) / 12;
    let julian = date.day() as i32 - 32075
        + 1461 * (date.year() + 4800 + l) / 4
        + 367 * (date.month() as i32 - 2 - l * 12) / 12
        - 3 * ((date.year() + 4900 + l) / 100) / 4;

    julian as f64 + (date.hour() as f64 / 24.0) - 0.5
}

/// Calculate days since j2000 using a UTC datetime object.
#[inline]
pub fn days_since_j2000(date: &DateTime<Utc>) -> f64 {
    const JD2000: f64 = 2451545.0;

    jdconv(date) - JD2000
}

#[inline]
pub fn spice_ephemeris_time(date: &DateTime<Utc>) -> f64 {
    let formatted_time = format!("{}", date.format("%Y %B %d, %H:%M:%S%.f %Z"));

    spice::str2et(formatted_time.as_str())
}

/// Convert civil time into local mean sidereal time.
/// Ref: https://idlastro.gsfc.nasa.gov/ftp/pro/astro/ct2lst.pro
///
/// ### Arguments
/// * 'lng' - The longitude in degrees (east of Greenwich) of the place for
///           which the local sidereal time is desired, scalar.   The Greenwich
///           mean sidereal time (GMST) can be found by setting Lng = 0.
/// * 'julian_day' -  Time of day of the specified date in decimal hours from j2000 epoch.
///
/// ### Returns
///     Mean sidereal time following input arguments as a double.
#[inline]
pub fn ct2lst(lng: f64, julian_day: f64) -> f64 {
    const C0: f64 = 280.46061837;
    const C1: f64 = 360.98564736629;
    const C2: f64 = 0.000387933;
    const C3: f64 = 38710000.0;
    const JULIAN_CENTURY: f64 = 36525.0;

    let t = julian_day / JULIAN_CENTURY;
    let theta = C0 + (C1 * julian_day) + t * t * (C2 - t / C3);

    ((theta + lng) / 15.0).rem_euclid(24.0)
}

#[cfg(test)]
mod tests {
    use chrono::TimeZone;

    use super::*;

    #[test]
    fn test_jdconv() {
        let jd_conversion = jdconv(&chrono::Utc.with_ymd_and_hms(1978, 1, 1, 0, 0, 0).unwrap());
        assert_eq!(jd_conversion, 2443509.5);
    }

    #[test]
    fn test_days_since_j2000() {
        let days = days_since_j2000(&chrono::Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap());
        assert_eq!(0.0, days);
    }

    #[test]
    fn test_ct2lst() {
        let lst = ct2lst(
            -76.72,
            days_since_j2000(
                &chrono::Utc
                    .with_ymd_and_hms(2008, 7, 30, 20, 53, 00)
                    .unwrap(),
            ),
        );

        assert!(11.356505 < lst && lst < 11.556505);
    }
}
