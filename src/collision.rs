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
    // Overlap checker for well formed ranges.
    let does_overlap = |a: (f64, f64), b: (f64, f64)| a.0.max(b.0) <= a.1.min(b.1);
    // This closure is used to fix improper ranges such that the property a.0 <= a.1 is maintained.
    let well_formed_range = |a: (f64, f64)| {
        if a.0 > a.1 {
            (a.1, a.0)
        } else {
            a
        }
    };
    for sim_body_i in 1..sim_bodies.len() {
        let back_range = {
            if let Some(range) = group_range {
                range
            } else {
                well_formed_range(accessor(sim_bodies.get(sim_body_i - 1).unwrap()))
            }
        };

        let front_range = well_formed_range(accessor(sim_bodies.get_mut(sim_body_i).unwrap()));
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

fn mark_overlapping_groups_axis<F, C>(
    sim_bodies: &mut [SimobjT],
    accessor: F,
    coord_cmp: C,
) -> usize
where
    F: Fn(&SimobjT) -> (f64, f64),
    C: Fn(&SimobjT, &SimobjT) -> std::cmp::Ordering + std::marker::Sync,
{
    let mut num_groups: usize = 0;

    for slice in sim_bodies.chunk_by_mut(|a, b| a.overlap_marker == b.overlap_marker) {
        // Skip over non-overlap group.
        if slice.first().unwrap().overlap_marker.is_none() {
            continue;
        }
        // Sort the slice by coord to enable building up of overlap groups. This allows for sequential
        // joining of simulation bodies into a given overlap group.
        slice.par_sort_unstable_by(&coord_cmp);
        // Unset previous overlap markers as to find new overlap groups for the provided axis. This is okay to
        // do here as the previous over lap groups has already been pulled out into a slice at this point.
        slice.iter_mut().for_each(|a| a.overlap_marker = None);
        num_groups = mark_overlapping_groups_axis_slice(slice, num_groups, &accessor);
    }

    // Return early if no overlap groups were found.
    if num_groups == 0 {
        return 0;
    }

    // Group up non-overlap set and overlapping groups for pulling out these segments in chunks.
    sim_bodies.par_sort_unstable_by(|a, b| a.overlap_marker.cmp(&b.overlap_marker));

    num_groups
}

pub fn find_collision_set(sim_bodies: &mut [SimobjT]) -> bool {
    // Place all objects within the same marker group initially (1) as the set has yet
    // to be reduced into discrete groups.
    sim_bodies
        .iter_mut()
        .for_each(|a| a.overlap_marker = Some(1));

    let x_accessor = |a: &SimobjT| (a.saved_state.coord_helio.x, a.state.coord_helio.x);
    let x_cmp = |a: &SimobjT, b: &SimobjT| {
        a.saved_state
            .coord_helio
            .x
            .partial_cmp(&b.saved_state.coord_helio.x)
            .unwrap()
    };
    if mark_overlapping_groups_axis(sim_bodies, x_accessor, x_cmp) == 0 {
        return false;
    }

    let y_accessor = |a: &SimobjT| (a.saved_state.coord_helio.y, a.state.coord_helio.y);
    let y_cmp = |a: &SimobjT, b: &SimobjT| {
        a.saved_state
            .coord_helio
            .y
            .partial_cmp(&b.saved_state.coord_helio.y)
            .unwrap()
    };
    if mark_overlapping_groups_axis(sim_bodies, y_accessor, y_cmp) == 0 {
        return false;
    }

    let z_accessor = |a: &SimobjT| (a.saved_state.coord_helio.z, a.state.coord_helio.z);
    let z_cmp = |a: &SimobjT, b: &SimobjT| {
        a.saved_state
            .coord_helio
            .z
            .partial_cmp(&b.saved_state.coord_helio.z)
            .unwrap()
    };
    if mark_overlapping_groups_axis(sim_bodies, z_accessor, z_cmp) == 0 {
        return false;
    }

    // At this point the vector of simulation bodies are chunked into a group which does
    // not overlap and groups of overlapping bounding boxes numbered 1,2,3,4...
    // These groups flat within the array sorted by numbering. Groups may be pulled out
    // by inspecting the numbering (or lack thereof) via the 'overlap_marker' member.
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

    use super::{
        find_collision_set, mark_overlapping_groups_axis, mark_overlapping_groups_axis_slice,
    };

    fn make_test_sim_bodies() -> Vec<SimobjT> {
        let mut sim_obj_template = SimobjT {
            id: 1,
            sim_object_type: SimObjectType::Spacecraft,
            perturb_store: None,
            name: "test_object_1".to_string(),
            drag_area: 1.0,
            drag_coeff: 1.0,
            mass: 1.0,
            state: SimObjTState {
                soi: Solarobj::Earth,
                coords: Array3d::default(),
                velocity: Array3d::default(),
                coord_helio: Array3d {
                    x: 100.0,
                    y: 200.0,
                    z: 300.0,
                },
                coord_helio_previous: Array3d::default(),
                coords_fixed: LLH::default(),
            },
            saved_state: SimObjTState {
                soi: Solarobj::Earth,
                coords: Array3d::default(),
                velocity: Array3d::default(),
                coord_helio: Array3d {
                    x: 50.0,
                    y: 100.0,
                    z: 150.0,
                },
                coord_helio_previous: Array3d::default(),
                coords_fixed: LLH::default(),
            },
            overlap_marker: None,
        };
        let mut sim_bodies = Vec::new();
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.name = "test_object_2".to_string();
        sim_obj_template.state.coord_helio.x = -160.0;
        sim_obj_template.saved_state.coord_helio.x = 90.0;
        sim_obj_template.state.coord_helio.y = 0.0;
        sim_obj_template.saved_state.coord_helio.y = 1.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.name = "test_object_3".to_string();
        sim_obj_template.state.coord_helio.x = 250.0;
        sim_obj_template.saved_state.coord_helio.x = 90.0;
        sim_obj_template.state.coord_helio.y = 300.0;
        sim_obj_template.saved_state.coord_helio.y = 150.0;
        sim_obj_template.state.coord_helio.z = 350.0;
        sim_obj_template.saved_state.coord_helio.z = 200.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.name = "test_object_4".to_string();
        sim_obj_template.state.coord_helio.x = 1000.0;
        sim_obj_template.saved_state.coord_helio.x = 1050.0;
        sim_obj_template.state.coord_helio.y = -1.0;
        sim_obj_template.saved_state.coord_helio.y = 1.1;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.name = "test_object_5".to_string();
        sim_obj_template.state.coord_helio.x = 10000.0;
        sim_obj_template.saved_state.coord_helio.x = 5050.0;
        sim_obj_template.state.coord_helio.y = -10.0;
        sim_obj_template.saved_state.coord_helio.y = 100.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.name = "test_object_6".to_string();
        sim_obj_template.state.coord_helio.x = 20000.0;
        sim_obj_template.saved_state.coord_helio.x = 7050.0;
        sim_obj_template.state.coord_helio.y = 160.0;
        sim_obj_template.saved_state.coord_helio.y = 350.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_bodies
    }

    #[test]
    fn test_mark_overlapping_groups_axis_slice() {
        let mut sim_bodies = make_test_sim_bodies();

        let group_cnt = mark_overlapping_groups_axis_slice(&mut sim_bodies, 0, &|a| {
            (a.saved_state.coord_helio.x, a.state.coord_helio.x)
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
    fn test_mark_overlapping_groups_axis() {
        let mut sim_bodies = make_test_sim_bodies();

        sim_bodies
            .iter_mut()
            .for_each(|a| a.overlap_marker = Some(1));

        let x_accessor = |a: &SimobjT| (a.saved_state.coord_helio.x, a.state.coord_helio.x);
        let x_cmp = |a: &SimobjT, b: &SimobjT| {
            a.saved_state
                .coord_helio
                .x
                .partial_cmp(&b.saved_state.coord_helio.x)
                .unwrap()
        };
        let group_cnt = mark_overlapping_groups_axis(&mut sim_bodies, x_accessor, x_cmp);
        assert_eq!(group_cnt, 2);
        assert!(sim_bodies.get(0).unwrap().overlap_marker.is_none());
        assert_eq!(sim_bodies.get(1).unwrap().overlap_marker.unwrap(), 1);
        assert_eq!(sim_bodies.get(2).unwrap().overlap_marker.unwrap(), 1);
        assert_eq!(sim_bodies.get(3).unwrap().overlap_marker.unwrap(), 1);
        assert_eq!(sim_bodies.get(4).unwrap().overlap_marker.unwrap(), 2);
        assert_eq!(sim_bodies.get(5).unwrap().overlap_marker.unwrap(), 2);

        let y_accessor = |a: &SimobjT| (a.saved_state.coord_helio.y, a.state.coord_helio.y);
        let y_cmp = |a: &SimobjT, b: &SimobjT| {
            a.saved_state
                .coord_helio
                .y
                .partial_cmp(&b.saved_state.coord_helio.y)
                .unwrap()
        };
        let group_cnt = mark_overlapping_groups_axis(&mut sim_bodies, y_accessor, y_cmp);
        assert_eq!(group_cnt, 1);
        assert!(sim_bodies.get(0).unwrap().overlap_marker.is_none());
        assert!(sim_bodies.get(1).unwrap().overlap_marker.is_none());
        assert!(sim_bodies.get(2).unwrap().overlap_marker.is_none());
        assert!(sim_bodies.get(3).unwrap().overlap_marker.is_none());
        assert_eq!(sim_bodies.get(4).unwrap().overlap_marker.unwrap(), 1);
        assert_eq!(sim_bodies.get(5).unwrap().overlap_marker.unwrap(), 1);
    }

    #[test]
    fn test_find_collision_set() {
        let mut sim_bodies = make_test_sim_bodies();
        let overlap_found = find_collision_set(&mut sim_bodies);

        assert!(overlap_found);
        assert!(sim_bodies.get(0).unwrap().overlap_marker.is_none());
        assert!(sim_bodies.get(1).unwrap().overlap_marker.is_none());
        assert!(sim_bodies.get(2).unwrap().overlap_marker.is_none());
        assert!(sim_bodies.get(3).unwrap().overlap_marker.is_none());
        assert_eq!(sim_bodies.get(4).unwrap().overlap_marker.unwrap(), 1);
        assert_eq!(sim_bodies.get(5).unwrap().overlap_marker.unwrap(), 1);
    }
}
