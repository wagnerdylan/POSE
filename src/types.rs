use std::ops;

#[derive(Debug, Clone, Copy)]
pub struct Array3d {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Array3d {
    pub fn dot(&self, other: &Array3d) -> f64 {
        (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
    }
}

impl ops::Add<Array3d> for Array3d {
    type Output = Array3d;

    fn add(self, rhs: Self) -> Self::Output {
        Array3d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl ops::Add<f64> for Array3d {
    type Output = Array3d;

    fn add(self, rhs: f64) -> Self::Output {
        Array3d {
            x: self.x + rhs,
            y: self.y + rhs,
            z: self.z + rhs,
        }
    }
}

impl ops::Mul<Array3d> for Array3d {
    type Output = Array3d;

    fn mul(self, rhs: Array3d) -> Self::Output {
        Array3d {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
    }
}

impl ops::Mul<f64> for Array3d {
    type Output = Array3d;

    fn mul(self, rhs: f64) -> Self::Output {
        Array3d {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}