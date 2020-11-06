#[derive(Debug, Clone, Copy)]
pub struct Array3d {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Array3d {
    pub fn add(&self, other: &Array3d) -> Array3d {
        Array3d {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    pub fn dot(&self, other: &Array3d) -> f64 {
        (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
    }
}
