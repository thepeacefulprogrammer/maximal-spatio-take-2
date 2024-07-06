import time as tm
from collections import defaultdict
from itertools import combinations
import bisect
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from algorithms.utils import manhattan_distance, chebyshev_distance


def rotate_point_pi_over_4(point):
    q1 = point.x
    q2 = point.y
    rotated_q1 = (q1 - q2) / math.sqrt(2)
    rotated_q2 = (q1 + q2) / math.sqrt(2)
    return rotated_q1, rotated_q2

def object_bounding_box_manhattan(position, d):
    q1, q2 = position
    min_point = (q1 - d, q2 - d)
    max_point = (q1 + d, q2 + d)
    return min_point, max_point

def object_bounding_box_chebyshev(position, d):
    q1, q2 = position
    min_point = (q1 - d, q2 - d)
    max_point = (q1 + d, q2 + d)
    return min_point, max_point

def bounding_box(points):
    min_point = tuple(map(min, zip(*points)))
    max_point = tuple(map(max, zip(*points)))
    return min_point, max_point

def composite_bounding_box(points, d, distance_fn):
    rotated_points = [rotate_point_pi_over_4(p) for p in points]
    bounding_boxes = [object_bounding_box_manhattan(p, d) if distance_fn == manhattan_distance else object_bounding_box_chebyshev(p, d) for p in rotated_points]
    min_points, max_points = zip(*bounding_boxes)
    return bounding_box(min_points), bounding_box(max_points)

class CompositeBoundingBox:
    def __init__(self, points, distance, distance_func):
        self.min_coords, self.max_coords = composite_bounding_box(points, distance, distance_func)

    def intersects(self, other):
        # Check AABB intersection
        return not (self.max_coords[0] < other.min_coords[0] or
                    self.min_coords[0] > other.max_coords[0] or
                    self.max_coords[1] < other.min_coords[1] or
                    self.min_coords[1] > other.max_coords[1])


def is_prevalent(candidate, instance_counts, bounding_boxes, minprev, intersection_cache):
    feature_counts = {feature: instance_counts[feature] for feature in candidate}
    if len(feature_counts) != len(candidate):
        return False

    participation_counts = defaultdict(int)
    for feature in candidate:
        sorted_instances = bounding_boxes[feature]
        other_features = candidate - {feature}

        for other_feature in other_features:
            if (feature, other_feature) in intersection_cache:
                intersecting_instances = intersection_cache[(feature, other_feature)]
            else:
                intersecting_instances = set()
                other_instances = bounding_boxes[other_feature]
                i, j = 0, 0
                while i < len(sorted_instances) and j < len(other_instances):
                    if sorted_instances[i].max_coords[0] < other_instances[j].min_coords[0]:
                        i += 1
                    elif other_instances[j].max_coords[0] < sorted_instances[i].min_coords[0]:
                        j += 1
                    else:
                        if sorted_instances[i].intersects(other_instances[j]):
                            intersecting_instances.add((i, j))
                        i += 1
                intersection_cache[(feature, other_feature)] = intersecting_instances

            participation_counts[feature] += len(intersecting_instances)
            if participation_counts[feature] >= feature_counts[feature] * minprev:
                break

    prevalences = [participation_counts[feature] / feature_counts[feature] for feature in candidate]
    return min(prevalences) >= minprev if prevalences else False


def binary_search(sorted_list, target, key=lambda x: x):
    return bisect.bisect_left(sorted_list, target, key=key)

def process_time_step_bounding_boxes(time, state, maxdist, distance_metric, all_features, bounding_boxes_per_time, cache):
    current_features = set(inst.id.feature for inst in state.instances)
    all_features.update(current_features)
    bounding_boxes = defaultdict(list)
    
    for feature in current_features:
        if feature in cache:
            bounding_boxes[feature] = cache[feature]
        else:
            feature_points = [inst.pos for inst in state.instances if inst.id.feature == feature]
            if feature_points:
                bb = CompositeBoundingBox(feature_points, maxdist, distance_metric)
                bounding_boxes[feature].append(bb)
            cache[feature] = bounding_boxes[feature]

    for feature in bounding_boxes:
        bounding_boxes[feature].sort(key=lambda bb: bb.min_coords[0])

    bounding_boxes_per_time[time] = bounding_boxes

def build_bounding_boxes(input_data, maxdist, distance_metric):
    all_features = set()
    bounding_boxes_per_time = {}
    cache = {}

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_time_step_bounding_boxes, time, state, maxdist, distance_metric, all_features, bounding_boxes_per_time, cache) for time, state in input_data.states.items()]
        for future in as_completed(futures):
            future.result()

    return all_features, bounding_boxes_per_time, cache

def process_time_step_prevalent_pairs(time, state, bounding_boxes_per_time, minprev, prevalent_pairs, time_prevalent_pairs, intersection_cache):
    bounding_boxes = bounding_boxes_per_time[time]
    instance_counts = state.count_instances()
    current_features = set(inst.id.feature for inst in state.instances)
    
    sorted_features = sorted(current_features)
    for i, f1 in enumerate(sorted_features):
        for f2 in sorted_features[i + 1:]:
            pair = frozenset([f1, f2])
            if is_prevalent(pair, instance_counts, bounding_boxes, minprev, intersection_cache):
                prevalent_pairs[time].add(pair)
                time_prevalent_pairs.add(pair)

def find_prevalent_pairs(all_features, bounding_boxes_per_time, input_data, minprev):
    prevalent_pairs = defaultdict(set)
    time_prevalent_pairs = set()
    intersection_cache = {}

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_time_step_prevalent_pairs, time, state, bounding_boxes_per_time, minprev, prevalent_pairs, time_prevalent_pairs, intersection_cache) for time, state in input_data.states.items()]
        for future in as_completed(futures):
            future.result()

    return prevalent_pairs, time_prevalent_pairs


def step1a_build_bounding_boxes(input_data, maxdist, distance_metric):
    return build_bounding_boxes(input_data, maxdist, distance_metric)

def step1b_find_prevalent_pairs(all_features, bounding_boxes_per_time, input_data, minprev):
    return find_prevalent_pairs(all_features, bounding_boxes_per_time, input_data, minprev)

def generate_candidates(prev_candidates, k):
    new_candidates = set()
    sorted_prev = sorted(prev_candidates)
    for i, c1 in enumerate(sorted_prev):
        for c2 in sorted_prev[i+1:]:
            new_candidate = c1.union(c2)
            if len(new_candidate) == k:
                if all(frozenset(comb) in prev_candidates for comb in combinations(new_candidate, k-1)):
                    new_candidates.add(new_candidate)
    return new_candidates

def step2_find_time_prevalent_pairs(prevalent_pairs, time_prevalent_pairs, minfreq, num_states):
    return {
        pair for pair in time_prevalent_pairs
        if sum(pair in pairs for pairs in prevalent_pairs.values()) >= minfreq * num_states
    }

def step3_generate_candidates(time_prevalent_pairs):
    candidates = list(time_prevalent_pairs)
    k = 3  # Start from size 3 candidates
    all_candidates = []

    while candidates:
        new_candidates = generate_candidates(candidates, k)
        all_candidates.append(new_candidates)
        candidates = new_candidates
        k += 1
    
    return all_candidates

def step4_prune_candidates(input_data, all_candidates, bounding_boxes_per_time, minfreq, minprev, prevalent_pairs):
    prevalence_cache = {}
    result = list(all_candidates[0])  # Start with prevalent pairs

    for candidates in all_candidates[1:]:
        prevalent_candidates = []

        with ThreadPoolExecutor() as executor:
            future_to_candidate = {executor.submit(check_prevalence, candidate, bounding_boxes_per_time, input_data, minfreq, minprev, prevalence_cache, prevalent_pairs): candidate for candidate in candidates}
            for future in as_completed(future_to_candidate):
                candidate = future_to_candidate[future]
                is_prevalent = future.result()
                if is_prevalent:
                    prevalent_candidates.append(candidate)

        result.extend(prevalent_candidates)
    
    return result

def check_prevalence(candidate, bounding_boxes_per_time, input_data, minfreq, minprev, cache, prevalent_pairs):
    if candidate in cache:
        return cache[candidate]

    time_prevalent_count = 0
    for time, bounding_boxes in bounding_boxes_per_time.items():
        state = input_data.states[time]
        instance_counts = state.count_instances()
        if all(feature in instance_counts for feature in candidate):
            if improved_instance_identification(candidate, bounding_boxes_per_time, state, time, minprev, cache):
                time_prevalent_count += 1
                if time_prevalent_count >= minfreq * len(input_data.states):
                    cache[candidate] = True
                    return True
    cache[candidate] = False
    return False

def improved_instance_identification(candidate, bounding_boxes_per_time, state, time, minprev, cache):
    if candidate in cache:
        return cache[candidate]
    
    instance_counts = state.count_instances()
    bounding_boxes = bounding_boxes_per_time[time]
    sorting_dimension = 0  # Assuming sorting by the first dimension for binary search

    for feature in candidate:
        if feature not in instance_counts:
            cache[candidate] = False
            return False
        
        other_features = candidate - {feature}
        sorted_instances = sorted(bounding_boxes[feature], key=lambda bb: bb.min_coords[sorting_dimension])

        for other_feature in other_features:
            count = 0
            for bb in bounding_boxes[other_feature]:
                index = binary_search(sorted_instances, bb.min_coords[sorting_dimension], key=lambda x: x.max_coords[sorting_dimension])
                if index < len(sorted_instances):
                    upper_bound = bb.max_coords[sorting_dimension]
                    for inst in sorted_instances[index:]:
                        if inst.min_coords[sorting_dimension] > upper_bound:
                            break
                        if inst.intersects(bb):
                            count += 1
                            if count / instance_counts[feature] >= minprev:
                                break
            if count / instance_counts[feature] < minprev:
                cache[candidate] = False
                return False
    cache[candidate] = True
    return True

def step5_ensure_maximality(result):
    maximal_result = []
    for pattern in sorted(result, key=len, reverse=True):
        if not any(pattern.issubset(maximal) for maximal in maximal_result):
            maximal_result.append(pattern)
    return maximal_result

def BBmaxspatiotempcolloc(input_data, maxdist, minprev, minfreq, distance_metric=chebyshev_distance, predictions=False, verbose=0, clean_trees=False):
    start_time = tm.time()

    # Step 1A
    step1a_start = tm.time()
    all_features, bounding_boxes_per_time, cache = step1a_build_bounding_boxes(input_data, maxdist, distance_metric)
    step1a_end = tm.time()
    if verbose > 0:
        print(f"Step 1A (Build bounding boxes) took {step1a_end - step1a_start:.2f} seconds")

    # Step 1B
    step1b_start = tm.time()
    prevalent_pairs, time_prevalent_pairs = step1b_find_prevalent_pairs(all_features, bounding_boxes_per_time, input_data, minprev)
    step1b_end = tm.time()
    if verbose > 0:
        print(f"Step 1B (Find prevalent pairs) took {step1b_end - step1b_start:.2f} seconds")

    # Step 2
    step2_start = tm.time()
    time_prevalent_pairs = step2_find_time_prevalent_pairs(prevalent_pairs, time_prevalent_pairs, minfreq, len(input_data.states))
    step2_end = tm.time()
    if verbose > 0:
        print(f"Step 2 (Find time-prevalent pairs) took {step2_end - step2_start:.2f} seconds")

    if not time_prevalent_pairs:
        return []

    # Step 3
    step3_start = tm.time()
    all_candidates = step3_generate_candidates(time_prevalent_pairs)
    step3_end = tm.time()
    if verbose > 0:
        print(f"Step 3 (Generate candidates) took {step3_end - step3_start:.2f} seconds")

    # Step 4
    step4_start = tm.time()
    result = step4_prune_candidates(input_data, all_candidates, bounding_boxes_per_time, minfreq, minprev, prevalent_pairs)
    step4_end = tm.time()
    if verbose > 0:
        print(f"Step 4 (Prune candidates) took {step4_end - step4_start:.2f} seconds")

    # Step 5
    step5_start = tm.time()
    maximal_result = step5_ensure_maximality(result)
    step5_end = tm.time()
    if verbose > 0:
        print(f"Step 5 (Ensure maximality) took {step5_end - step5_start:.2f} seconds")

    end_time = tm.time()
    if verbose > 0:
        print(f"Total execution time: {end_time - start_time:.2f} seconds")
    
    return maximal_result
