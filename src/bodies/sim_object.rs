use super::solar_model::Solarobj;
use crate::{
    output::{self, SimulationObjectPerturbationOut},
    types::{Array3d, LLH},
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Default)]
pub enum SimObjectType {
    #[default]
    Spacecraft,
    Debris,
}

pub struct PerturbationDefinition {
    pub perturb_name: String,
    pub x_accel: f64,
    pub y_accel: f64,
    pub z_accel: f64,
}

#[derive(Default)]
pub struct PerturbationStore {
    pub sun_gravity: Option<PerturbationDefinition>,
    pub earth_gravity: Option<PerturbationDefinition>,
    pub moon_gravity: Option<PerturbationDefinition>,
    pub earth_drag: Option<PerturbationDefinition>,
}

impl PerturbationStore {
    fn make_output_object(
        id: u32,
        name: &str,
        sim_time: f64,
        definition: &PerturbationDefinition,
    ) -> SimulationObjectPerturbationOut {
        SimulationObjectPerturbationOut {
            id,
            sim_time,
            name: name.to_string(),
            perturb_type: definition.perturb_name.clone(),
            x_accel: definition.x_accel,
            y_accel: definition.y_accel,
            z_accel: definition.z_accel,
        }
    }

    pub fn to_output_form(
        &self,
        id: u32,
        name: &str,
        sim_time: f64,
    ) -> Vec<SimulationObjectPerturbationOut> {
        let mut out = Vec::new();

        if let Some(perturb) = &self.sun_gravity {
            out.push(PerturbationStore::make_output_object(
                id, name, sim_time, perturb,
            ))
        }

        if let Some(perturb) = &self.earth_gravity {
            out.push(PerturbationStore::make_output_object(
                id, name, sim_time, perturb,
            ))
        }

        if let Some(perturb) = &self.moon_gravity {
            out.push(PerturbationStore::make_output_object(
                id, name, sim_time, perturb,
            ))
        }

        if let Some(perturb) = &self.earth_drag {
            out.push(PerturbationStore::make_output_object(
                id, name, sim_time, perturb,
            ))
        }

        out
    }
}

#[derive(Serialize, Deserialize)]
pub struct SimobjT {
    #[serde(skip_deserializing)]
    pub id: u32,
    #[serde(skip_deserializing)]
    pub sim_object_type: SimObjectType,
    #[serde(skip_deserializing)]
    pub coords_abs: Array3d,
    #[serde(skip_deserializing)]
    pub coords_fixed: LLH,
    #[serde(skip_deserializing, skip_serializing)]
    pub perturb_store: Option<PerturbationStore>,
    pub name: String,
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
        let fixed_coords = &self.coords_fixed;
        let coords = &self.coords;
        let velocity = &self.velocity;

        output::SimulationObjectParameters {
            id: self.id,
            sim_time,
            name: self.name.clone(),
            soi: self.soi.to_string(),
            x_abs_coord: abs_coords.x,
            y_abs_coord: abs_coords.y,
            z_abs_coord: abs_coords.z,
            lat: fixed_coords.lat,
            long: fixed_coords.long,
            altitude: fixed_coords.alt,
            x_coord: coords.x,
            y_coord: coords.y,
            z_coord: coords.z,
            x_velocity: velocity.x,
            y_velocity: velocity.y,
            z_velocity: velocity.z,
        }
    }
}
