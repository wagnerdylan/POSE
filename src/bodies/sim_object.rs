use super::solar_object::Solarobj;
use crate::{output, types::Array3d};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Default)]
pub enum SimObjectType {
    #[default]
    Spacecraft,
    Debris,
}

#[derive(Serialize, Deserialize)]
pub struct SimobjT {
    #[serde(skip_deserializing)]
    pub id: u32,
    #[serde(skip_deserializing)]
    pub sim_object_type: SimObjectType,
    #[serde(skip_deserializing)]
    pub coords_abs: Array3d,
    pub soi: Solarobj,
    pub coords: Array3d,
    pub velocity: Array3d,
    pub drag_area: f64,
    pub drag_coeff: f64,
    pub mass: f64,
}

impl SimobjT {
    pub fn to_output_form(&self, sim_time: f64) -> output::SimulationObjectParameters {
        let abs_coords = &self.coords_abs;
        let coords = &self.coords;
        let velocity = &self.velocity;

        output::SimulationObjectParameters {
            id: self.id,
            sim_time,
            soi: self.soi.to_string(),
            x_abs_coord: abs_coords.x,
            y_abs_coord: abs_coords.y,
            z_abs_coord: abs_coords.z,
            x_coord: coords.x,
            y_coord: coords.y,
            z_coord: coords.z,
            x_velocity: velocity.x,
            y_velocity: velocity.y,
            z_velocity: velocity.z,
        }
    }
}
