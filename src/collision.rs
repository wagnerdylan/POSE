use rayon::slice::ParallelSliceMut;

use crate::bodies::sim_object::SimobjT;

fn mark_overlapping_groups_axis_slice<F>(
    sim_bodies: &mut [SimobjT],
    mut group_cnt: usize,
    accessor: &F,
) -> usize
where
    F: Fn(&SimobjT) -> (f64, f64),
{
    let mut group_range: Option<(f64, f64)> = None;
    let does_overlap = |a: (f64, f64), b: (f64, f64)| a.0.max(b.0) <= a.1.min(a.1);
    for sim_body_i in 1..sim_bodies.len() {
        let back_range = {
            if let Some(range) = group_range {
                range
            } else {
                accessor(sim_bodies.get(sim_body_i - 1).unwrap())
            }
        };

        let front_range = accessor(sim_bodies.get_mut(sim_body_i).unwrap());
        if does_overlap(back_range, front_range) {
            let marker_num = {
                if let Some(marker) = sim_bodies.get_mut(sim_body_i - 1).unwrap().overlap_marker {
                    marker
                } else {
                    group_cnt += 1;
                    group_cnt
                }
            };
            sim_bodies.get_mut(sim_body_i - 1).unwrap().overlap_marker = Some(marker_num);
            sim_bodies.get_mut(sim_body_i).unwrap().overlap_marker = Some(marker_num);
            group_range = Some((back_range.0, front_range.1))
        } else {
            group_range = None;
        }
    }

    group_cnt
}

fn mark_overlapping_groups_axis<F>(sim_bodies: &mut [SimobjT], accessor: F) -> usize
where
    F: Fn(&SimobjT) -> (f64, f64),
{
    let mut num_groups: usize = 0;

    for slice in sim_bodies.chunk_by_mut(|a, b| a.overlap_marker == b.overlap_marker) {
        // Skip over non-collision group.
        if slice.first().unwrap().overlap_marker.is_none() {
            continue;
        }
        slice.iter_mut().for_each(|a| a.overlap_marker = None);
        num_groups = mark_overlapping_groups_axis_slice(slice, num_groups, &accessor);
    }

    // Return early if no collision groups were found.
    if num_groups == 0 {
        return 0;
    }

    sim_bodies.par_sort_unstable_by(|a, b| a.overlap_marker.cmp(&b.overlap_marker));

    num_groups
}

pub fn find_collision_set(sim_bodies: &mut [SimobjT]) -> bool {
    // Sort sim_bodies list by the x axis to prepare for overlapping groups search.
    sim_bodies.par_sort_unstable_by(|a, b| {
        a.saved_state
            .coords_abs
            .x
            .partial_cmp(&b.saved_state.coords_abs.x)
            .unwrap()
    });
    // Place all objects within the same marker group initially (1) as the set has yet
    // to be reduced into discrete groups.
    sim_bodies
        .iter_mut()
        .for_each(|a| a.overlap_marker = Some(1));
    if mark_overlapping_groups_axis(sim_bodies, |a| {
        (a.saved_state.coords_abs.x, a.state.coords_abs.x)
    }) == 0
    {
        return false;
    }
    if mark_overlapping_groups_axis(sim_bodies, |a| {
        (a.saved_state.coords_abs.y, a.state.coords_abs.y)
    }) == 0
    {
        return false;
    }
    if mark_overlapping_groups_axis(sim_bodies, |a| {
        (a.saved_state.coords_abs.z, a.state.coords_abs.z)
    }) == 0
    {
        return false;
    }

    true
}

#[cfg(test)]
mod tests {
    use crate::{
        bodies::{
            sim_object::{SimObjTState, SimObjectType, SimobjT},
            solar_model::Solarobj,
        },
        types::{Array3d, LLH},
    };

    use super::mark_overlapping_groups_axis_slice;

    #[test]
    fn test_mark_overlapping_groups_axis_slice() {
        let mut sim_obj_template = SimobjT {
            id: 1,
            sim_object_type: SimObjectType::Spacecraft,
            perturb_store: None,
            name: "test_object".to_string(),
            drag_area: 1.0,
            drag_coeff: 1.0,
            mass: 1.0,
            state: SimObjTState {
                soi: Solarobj::Earth,
                coords: Array3d::default(),
                velocity: Array3d::default(),
                coords_abs: Array3d {
                    x: 100.0,
                    y: 200.0,
                    z: 300.0,
                },
                coords_abs_previous: Array3d::default(),
                coords_fixed: LLH::default(),
            },
            saved_state: SimObjTState {
                soi: Solarobj::Earth,
                coords: Array3d::default(),
                velocity: Array3d::default(),
                coords_abs: Array3d {
                    x: 50.0,
                    y: 100.0,
                    z: 150.0,
                },
                coords_abs_previous: Array3d::default(),
                coords_fixed: LLH::default(),
            },
            overlap_marker: None,
        };
        let mut sim_bodies = Vec::new();
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.state.coords_abs.x = 160.0;
        sim_obj_template.saved_state.coords_abs.x = 100.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.state.coords_abs.x = 250.0;
        sim_obj_template.saved_state.coords_abs.x = 90.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.state.coords_abs.x = 1000.0;
        sim_obj_template.saved_state.coords_abs.x = 1050.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.state.coords_abs.x = 10000.0;
        sim_obj_template.saved_state.coords_abs.x = 5050.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.state.coords_abs.x = 20000.0;
        sim_obj_template.saved_state.coords_abs.x = 7050.0;
        sim_bodies.push(sim_obj_template.clone());

        let group_cnt = mark_overlapping_groups_axis_slice(&mut sim_bodies, 0, &|a| {
            (a.saved_state.coords_abs.x, a.state.coords_abs.x)
        });

        assert_eq!(group_cnt, 2);
        assert_eq!(sim_bodies.get(0).unwrap().overlap_marker.unwrap(), 1);
        assert_eq!(sim_bodies.get(1).unwrap().overlap_marker.unwrap(), 1);
        assert_eq!(sim_bodies.get(2).unwrap().overlap_marker.unwrap(), 1);
        assert!(sim_bodies.get(3).unwrap().overlap_marker.is_none());
        assert_eq!(sim_bodies.get(4).unwrap().overlap_marker.unwrap(), 2);
        assert_eq!(sim_bodies.get(5).unwrap().overlap_marker.unwrap(), 2);
    }

    #[test]
    fn test_mark_overlapping_groups_axis() {}

    #[test]
    fn find_collision_set() {}
}
