use super::solar_model::Solarobj;
use crate::{
    output::{self, SimulationObjectPerturbationOut},
    types::{Array3d, LLH},
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Default, Clone)]
pub enum SimObjectType {
    #[default]
    Spacecraft,
    Debris,
}

#[derive(Clone)]
pub struct PerturbationDefinition {
    pub perturb_name: String,
    pub x_accel: f64,
    pub y_accel: f64,
    pub z_accel: f64,
}

#[derive(Default, Clone)]
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

#[derive(Serialize, Deserialize, Default, Clone)]
pub struct SimObjTState {
    pub soi: Solarobj,
    pub coords: Array3d,
    pub velocity: Array3d,
    #[serde(skip_deserializing)]
    pub coord_helio: Array3d,
    #[serde(skip_deserializing)]
    pub coord_helio_previous: Array3d,
    #[serde(skip_deserializing)]
    pub coords_fixed: LLH,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct SimobjT {
    #[serde(skip_deserializing)]
    pub id: u32,
    #[serde(skip_deserializing)]
    pub sim_object_type: SimObjectType,
    #[serde(skip_deserializing, skip_serializing)]
    pub perturb_store: Option<PerturbationStore>,
    pub name: String,
    pub drag_area: f64,
    pub drag_coeff: f64,
    pub mass: f64,
    pub state: SimObjTState,
    #[serde(skip_deserializing, skip_serializing)]
    pub saved_state: SimObjTState,
    #[serde(skip_deserializing, skip_serializing)]
    pub overlap_marker: Option<usize>,
}

impl SimobjT {
    pub fn to_output_form(&self, sim_time: f64) -> output::SimulationObjectParameters {
        let coord_helio = &self.state.coord_helio;
        let fixed_coords = &self.state.coords_fixed;
        let coords = &self.state.coords;
        let velocity = &self.state.velocity;

        output::SimulationObjectParameters {
            id: self.id,
            sim_time,
            name: self.name.clone(),
            soi: self.state.soi.to_string(),
            x_coord_helio: coord_helio.x,
            y_coord_helio: coord_helio.y,
            z_coord_helio: coord_helio.z,
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
