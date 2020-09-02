pub struct Array3d {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Array3d {
    fn add(&self, other: &Array3d) -> Array3d {
        Array3d {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}
