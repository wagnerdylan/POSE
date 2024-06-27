use std::collections::HashSet;

use rayon::slice::ParallelSliceMut;

use crate::{
    bodies::sim_object::{SimObjectType, SimobjT},
    environment::Environment,
    output::{self, CollisionInfoOut},
    types::{self, Array3d},
};

pub struct IntersectionResult {
    pub body_a_idx: usize,
    pub body_b_idx: usize,
    pub body_a_intersect_coord: Array3d,
    pub body_b_intersect_coord: Array3d,
    pub intersect_dist: f64,
}

pub struct CollisionResult {
    pub body_a_id: u32,
    pub body_b_id: u32,
    pub body_a_name: String,
    pub body_b_name: String,
    pub body_a_coord: Array3d,
    pub body_b_coord: Array3d,
    pub sim_time: f64,
    pub intercept_distance: f64,
    pub relative_velocity: f64,
    pub new_sim_bodies: Vec<SimobjT>,
}

impl CollisionResult {
    pub fn to_output_form(&self, new_bodies: &[SimobjT]) -> output::CollisionInfoOut {
        let mut bodies_ids = Vec::new();
        for body in new_bodies {
            bodies_ids.push(body.id.to_string());
        }
        let bodies_id_str = bodies_ids.join(" ");

        CollisionInfoOut {
            body_a_id: self.body_a_id,
            body_b_id: self.body_b_id,
            body_a_name: self.body_a_name.clone(),
            body_b_name: self.body_b_name.clone(),
            sim_time: self.sim_time,
            intercept_distance: self.intercept_distance,
            relative_velocity: self.relative_velocity,
            body_a_x_coord: self.body_a_coord.x,
            body_a_y_coord: self.body_a_coord.y,
            body_a_z_coord: self.body_a_coord.z,
            body_b_x_coord: self.body_b_coord.x,
            body_b_y_coord: self.body_b_coord.y,
            body_b_z_coord: self.body_b_coord.z,
            generated_bodies_id: bodies_id_str,
        }
    }
}

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

    // Group overlapping objects by the overlap marker to collect previously marked together.
    // Additionally group by SOI to handle edge cases where objects span multiple SOIs.
    for slice in sim_bodies.chunk_by_mut(|a, b| {
        a.overlap_marker == b.overlap_marker
            && a.state.soi == b.state.soi
            && a.saved_state.soi == b.saved_state.soi
    }) {
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

fn keep_only_satellite_intersection_groups(sim_bodies: &mut [SimobjT]) -> bool {
    let mut sat_group_nums = HashSet::new();
    for slice in sim_bodies.chunk_by_mut(|a, b| a.overlap_marker == b.overlap_marker) {
        // Skip over non-overlap group.
        if slice.first().unwrap().overlap_marker.is_none() {
            continue;
        }

        for object in slice {
            if let SimObjectType::Spacecraft = object.sim_object_type {
                // Unwrap is okay here as only populated overlap markers are allowed at this point.
                sat_group_nums.insert(object.overlap_marker.unwrap());
                break;
            }
        }
    }

    for slice in sim_bodies.chunk_by_mut(|a, b| a.overlap_marker == b.overlap_marker) {
        // Skip over non-overlap group.
        if slice.first().unwrap().overlap_marker.is_none() {
            continue;
        }

        if sat_group_nums.contains(&slice.first().unwrap().overlap_marker.unwrap()) {
            continue;
        }

        slice.iter_mut().for_each(|x| x.overlap_marker = None);
    }

    // Resort sim_bodies to move newly de-marked overlap groups into the "None" group.
    sim_bodies.par_sort_unstable_by(|a, b| a.overlap_marker.cmp(&b.overlap_marker));

    !sat_group_nums.is_empty()
}

pub fn find_collision_set(sim_bodies: &mut [SimobjT], check_only_satellites: bool) -> bool {
    // Place all objects within the same marker group initially (1) as the set has yet
    // to be reduced into discrete groups.
    sim_bodies
        .iter_mut()
        .for_each(|a| a.overlap_marker = Some(1));

    let x_accessor = |a: &SimobjT| (a.saved_state.coords.x, a.state.coords.x);
    let x_cmp = |a: &SimobjT, b: &SimobjT| {
        a.saved_state
            .coords
            .x
            .partial_cmp(&b.saved_state.coords.x)
            .unwrap()
    };
    if mark_overlapping_groups_axis(sim_bodies, x_accessor, x_cmp) == 0 {
        return false;
    }

    let y_accessor = |a: &SimobjT| (a.saved_state.coords.y, a.state.coords.y);
    let y_cmp = |a: &SimobjT, b: &SimobjT| {
        a.saved_state
            .coords
            .y
            .partial_cmp(&b.saved_state.coords.y)
            .unwrap()
    };
    if mark_overlapping_groups_axis(sim_bodies, y_accessor, y_cmp) == 0 {
        return false;
    }

    let z_accessor = |a: &SimobjT| (a.saved_state.coords.z, a.state.coords.z);
    let z_cmp = |a: &SimobjT, b: &SimobjT| {
        a.saved_state
            .coords
            .z
            .partial_cmp(&b.saved_state.coords.z)
            .unwrap()
    };
    if mark_overlapping_groups_axis(sim_bodies, z_accessor, z_cmp) == 0 {
        return false;
    }

    if check_only_satellites && !keep_only_satellite_intersection_groups(sim_bodies) {
        return false;
    }

    // At this point the vector of simulation bodies are chunked into a group which does
    // not overlap and groups of overlapping bounding boxes numbered 1,2,3,4...
    // These groups flat within the array sorted by numbering. Groups may be pulled out
    // by inspecting the numbering (or lack thereof) via the 'overlap_marker' member.
    true
}

fn line_line_n_point_dist(
    l1: (&Array3d, &Array3d),
    l2: (&Array3d, &Array3d),
) -> (f64, Array3d, Array3d) {
    let mut min_dist = f64::INFINITY;
    let mut min_interp_point = f64::INFINITY;
    let mut min_p12 = Array3d::default();
    let mut min_p34 = Array3d::default();

    let interp_num = 6;
    for i in 0..=interp_num {
        let interp_point = (1.0 / interp_num as f64) * i as f64;
        let p12 = ((1f64 - interp_point) * l1.0) + interp_point * l1.1;
        let p34 = ((1f64 - interp_point) * l2.0) + interp_point * l2.1;

        let dist = types::l2_norm(&(p34 - p12));
        if dist < min_dist {
            min_dist = dist;
            min_interp_point = interp_point;
            min_p12 = p12;
            min_p34 = p34;
        }
    }

    if types::l2_norm(&(l1.0 - l1.1)) < 1.0 || types::l2_norm(&(l2.0 - l2.1)) < 1.0 {
        return (min_dist, min_p12, min_p34);
    }

    let p12_min = ((1f64 - min_interp_point) * l1.0) + min_interp_point * l1.1;
    let p34_min = ((1f64 - min_interp_point) * l2.0) + min_interp_point * l2.1;
    let mut left_dist_min = (f64::INFINITY, Array3d::default(), Array3d::default());
    let mut right_dist_min = (f64::INFINITY, Array3d::default(), Array3d::default());

    if min_interp_point > 0.0 + f64::EPSILON {
        let left_interp_point = min_interp_point - (1.0 / interp_num as f64);
        let p12_left = ((1f64 - left_interp_point) * l1.0) + left_interp_point * l1.1;
        let p34_left = ((1f64 - left_interp_point) * l2.0) + left_interp_point * l2.1;
        left_dist_min = line_line_n_point_dist((&p12_left, &p12_min), (&p34_left, &p34_min));
    }

    if min_interp_point < 1.0 - f64::EPSILON {
        let right_interp_point = min_interp_point + (1.0 / interp_num as f64);
        let p12_right = ((1f64 - right_interp_point) * l1.0) + right_interp_point * l1.1;
        let p34_right = ((1f64 - right_interp_point) * l2.0) + right_interp_point * l2.1;
        right_dist_min = line_line_n_point_dist((&p12_min, &p12_right), (&p34_min, &p34_right));
    }

    if min_dist < left_dist_min.0 && min_dist < right_dist_min.0 {
        (min_dist, min_p12, min_p34)
    } else if left_dist_min.0 < min_dist && left_dist_min.0 < right_dist_min.0 {
        left_dist_min
    } else {
        right_dist_min
    }
}

pub fn find_body_intersections(
    sim_bodies: &[SimobjT],
    current_step_count: u64,
) -> Vec<IntersectionResult> {
    let mut result_vec = Vec::new();
    for i in 0..sim_bodies.len() {
        for j in i + 1..sim_bodies.len() {
            // Skip over simulation bodies which have been deleted on current step.
            let body_a = sim_bodies.get(i).unwrap();
            let body_b = sim_bodies.get(j).unwrap();
            if let Some(body_stop) = body_a.marked_for_deletion_on_step {
                if body_stop <= current_step_count {
                    continue;
                }
            }
            if let Some(body_stop) = body_b.marked_for_deletion_on_step {
                if body_stop <= current_step_count {
                    continue;
                }
            }

            // Skipping over simulation bodies which are not in the same soi
            // prevents dis-jointed calculations when using velocity between the
            // two bodies. This may result in a false-negative collision check although
            // this case should be expectantly rare.
            if body_a.state.soi != body_b.state.soi {
                continue;
            }

            // Calculate the shortest line between two trajectory lines of body_a and body_b.
            let intersect_line = line_line_n_point_dist(
                (&body_a.state.previous_coords, &body_a.state.coords),
                (&body_b.state.previous_coords, &body_b.state.coords),
            );
            let intersect_dist = intersect_line.0;
            let body_a_intersect_coord = intersect_line.1;
            let body_b_intersect_coord = intersect_line.2;
            let intersect_bound = body_a.radius + body_b.radius;
            if intersect_dist < intersect_bound {
                result_vec.push(IntersectionResult {
                    body_a_idx: i,
                    body_b_idx: j,
                    body_a_intersect_coord,
                    body_b_intersect_coord,
                    intersect_dist,
                });
            }
        }
    }
    result_vec
}

pub fn collision_model(
    env: &Environment,
    sim_bodies: &mut [SimobjT],
    intersect_info: &IntersectionResult,
) -> CollisionResult {
    // TODO flush out this section with an empirical collision model.
    sim_bodies
        .get_mut(intersect_info.body_a_idx)
        .unwrap()
        .marked_for_deletion_on_step = Some(env.step_count);
    sim_bodies
        .get_mut(intersect_info.body_b_idx)
        .unwrap()
        .marked_for_deletion_on_step = Some(env.step_count);

    let body_a = sim_bodies.get(intersect_info.body_a_idx).unwrap();
    let body_b = sim_bodies.get(intersect_info.body_b_idx).unwrap();
    CollisionResult {
        new_sim_bodies: Vec::new(),
        body_a_id: body_a.id,
        body_b_id: body_b.id,
        body_a_name: body_a.name.clone(),
        body_b_name: body_b.name.clone(),
        body_a_coord: intersect_info.body_a_intersect_coord,
        body_b_coord: intersect_info.body_b_intersect_coord,
        sim_time: env.sim_time_s,
        intercept_distance: intersect_info.intersect_dist,
        relative_velocity: types::l2_norm(&(body_a.state.velocity - body_b.state.velocity)).abs(),
    }
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
            radius: 1.0,
            state: SimObjTState {
                soi: Solarobj::Earth,
                coords: Array3d {
                    x: 100.0,
                    y: 200.0,
                    z: 300.0,
                },
                previous_coords: Array3d::default(),
                velocity: Array3d::default(),
                coord_helio: Array3d::default(),
                coords_fixed: LLH::default(),
            },
            saved_state: SimObjTState {
                soi: Solarobj::Earth,
                velocity: Array3d::default(),
                coords: Array3d {
                    x: 50.0,
                    y: 100.0,
                    z: 150.0,
                },
                coords_fixed: LLH::default(),
                previous_coords: Array3d::default(),
                coord_helio: Array3d::default(),
            },
            overlap_marker: None,
            marked_for_deletion_on_step: None,
        };
        let mut sim_bodies = Vec::new();
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.name = "test_object_2".to_string();
        sim_obj_template.state.coords.x = -160.0;
        sim_obj_template.saved_state.coords.x = 90.0;
        sim_obj_template.state.coords.y = 0.0;
        sim_obj_template.saved_state.coords.y = 1.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.name = "test_object_3".to_string();
        sim_obj_template.state.coords.x = 250.0;
        sim_obj_template.saved_state.coords.x = 90.0;
        sim_obj_template.state.coords.y = 300.0;
        sim_obj_template.saved_state.coords.y = 150.0;
        sim_obj_template.state.coords.z = 350.0;
        sim_obj_template.saved_state.coords.z = 200.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.name = "test_object_4".to_string();
        sim_obj_template.state.coords.x = 1000.0;
        sim_obj_template.saved_state.coords.x = 1050.0;
        sim_obj_template.state.coords.y = -1.0;
        sim_obj_template.saved_state.coords.y = 1.1;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.name = "test_object_5".to_string();
        sim_obj_template.state.coords.x = 10000.0;
        sim_obj_template.saved_state.coords.x = 5050.0;
        sim_obj_template.state.coords.y = -10.0;
        sim_obj_template.saved_state.coords.y = 100.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_obj_template.name = "test_object_6".to_string();
        sim_obj_template.state.coords.x = 20000.0;
        sim_obj_template.saved_state.coords.x = 7050.0;
        sim_obj_template.state.coords.y = 160.0;
        sim_obj_template.saved_state.coords.y = 350.0;
        sim_bodies.push(sim_obj_template.clone());

        sim_bodies
    }

    #[test]
    fn test_mark_overlapping_groups_axis_slice() {
        let mut sim_bodies = make_test_sim_bodies();

        let group_cnt = mark_overlapping_groups_axis_slice(&mut sim_bodies, 0, &|a| {
            (a.saved_state.coords.x, a.state.coords.x)
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

        let x_accessor = |a: &SimobjT| (a.saved_state.coords.x, a.state.coords.x);
        let x_cmp = |a: &SimobjT, b: &SimobjT| {
            a.saved_state
                .coords
                .x
                .partial_cmp(&b.saved_state.coords.x)
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

        let y_accessor = |a: &SimobjT| (a.saved_state.coords.y, a.state.coords.y);
        let y_cmp = |a: &SimobjT, b: &SimobjT| {
            a.saved_state
                .coords
                .y
                .partial_cmp(&b.saved_state.coords.y)
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
        let overlap_found = find_collision_set(&mut sim_bodies, false);

        assert!(overlap_found);
        assert!(sim_bodies.get(0).unwrap().overlap_marker.is_none());
        assert!(sim_bodies.get(1).unwrap().overlap_marker.is_none());
        assert!(sim_bodies.get(2).unwrap().overlap_marker.is_none());
        assert!(sim_bodies.get(3).unwrap().overlap_marker.is_none());
        assert_eq!(sim_bodies.get(4).unwrap().overlap_marker.unwrap(), 1);
        assert_eq!(sim_bodies.get(5).unwrap().overlap_marker.unwrap(), 1);
    }
}
