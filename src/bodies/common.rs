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
