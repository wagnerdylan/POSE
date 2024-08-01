use serde::{Deserialize, Serialize};
use std::{iter::Sum, ops};

#[derive(Serialize, Deserialize, Debug, Clone, Copy, Default)]
pub struct Array3d {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Serialize, Deserialize, Default, Clone)]
pub struct LLH {
    pub lat: f64,
    pub long: f64,
    pub alt: f64,
}

impl Array3d {
    pub fn dot(&self, other: &Array3d) -> f64 {
        (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
    }

    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

pub fn l2_norm(x: &Array3d) -> f64 {
    x.dot(x).sqrt()
}

impl<'a> Sum<&'a Self> for Array3d {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
    {
        iter.fold(Self::default(), |acc, x| acc + x)
    }
}

impl_op_ex!(+ |a: &Array3d, b: &Array3d| -> Array3d {
    Array3d {
        x: a.x + b.x,
        y: a.y + b.y,
        z: a.z + b.z,
    }
});

impl_op_ex_commutative!(+ |a: &Array3d, b: &f64| -> Array3d {
    Array3d {
        x: a.x + b,
        y: a.y + b,
        z: a.z + b,
    }
});

impl_op_ex!(-|a: &Array3d, b: &Array3d| -> Array3d {
    Array3d {
        x: a.x - b.x,
        y: a.y - b.y,
        z: a.z - b.z,
    }
});

impl_op_ex!(-|a: &Array3d, b: &f64| -> Array3d {
    Array3d {
        x: a.x - b,
        y: a.y - b,
        z: a.z - b,
    }
});

impl_op_ex!(*|a: &Array3d, b: &Array3d| -> Array3d {
    Array3d {
        x: a.x * b.x,
        y: a.y * b.y,
        z: a.z * b.z,
    }
});

impl_op_ex_commutative!(*|a: &Array3d, b: &f64| -> Array3d {
    Array3d {
        x: a.x * b,
        y: a.y * b,
        z: a.z * b,
    }
});

impl_op_ex!(/|a: &Array3d, b: &Array3d| -> Array3d {
    Array3d {
        x: a.x / b.x,
        y: a.y / b.y,
        z: a.z / b.z,
    }
});

impl_op_ex!(/|a: &Array3d, b: &f64| -> Array3d {
    Array3d {
        x: a.x / b,
        y: a.y / b,
        z: a.z / b,
    }
});
