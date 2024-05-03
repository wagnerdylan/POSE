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
    let mut behind_range: Option<(f64, f64)> = None;
    let does_overlap = |a: (f64, f64), b: (f64, f64)| a.0.max(b.0) <= a.1.min(a.1);
    for sim_body_i in 1..sim_bodies.len() - 1 {
        let back_range = {
            if let Some(range) = Some(behind_range) {
                range.unwrap()
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
            behind_range = Some((back_range.0, front_range.1))
        } else {
            behind_range = None;
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
