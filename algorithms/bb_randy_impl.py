import time as tm
from collections import defaultdict
from itertools import combinations
import bisect
from algorithms.utils import manhattan_distance, chebyshev_distance
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import math


class RotatedBoundingBox:
    def __init__(self, min_coords, max_coords, angle):
        self.min_coords = min_coords
        self.max_coords = max_coords
        self.angle = angle
        self.aabb_min, self.aabb_max = self.calculate_aabb()

    def calculate_aabb(self):
        corners = [
            rotate_point(self.min_coords[0], self.min_coords[1], self.angle),
            rotate_point(self.min_coords[0], self.max_coords[1], self.angle),
            rotate_point(self.max_coords[0], self.min_coords[1], self.angle),
            rotate_point(self.max_coords[0], self.max_coords[1], self.angle)
        ]
        min_x = min(corner[0] for corner in corners)
        min_y = min(corner[1] for corner in corners)
        max_x = max(corner[0] for corner in corners)
        max_y = max(corner[1] for corner in corners)
        return (min_x, min_y), (max_x, max_y)

    def intersects(self, other):
        # Check AABB intersection
        return not (self.aabb_max[0] < other.aabb_min[0] or
                    self.aabb_min[0] > other.aabb_max[0] or
                    self.aabb_max[1] < other.aabb_min[1] or
                    self.aabb_min[1] > other.aabb_max[1])

def rotate_point(x, y, angle):
    rad = math.radians(angle)
    cos_rad = math.cos(rad)
    sin_rad = math.sin(rad)
    x_new = x * cos_rad - y * sin_rad
    y_new = x * sin_rad + y * cos_rad
    return x_new, y_new

def create_rotated_bounding_box(instance, distance, distance_func, angle):
    if distance_func == chebyshev_distance:
        offset = distance
    elif distance_func == manhattan_distance:
        offset = distance / 2
    else:
        raise ValueError("Unsupported distance function")

    min_coords = [instance.pos.x - offset, instance.pos.y - offset]
    max_coords = [instance.pos.x + offset, instance.pos.y + offset]

    return RotatedBoundingBox(min_coords, max_coords, angle)

def binary_search(sorted_list, target, key=lambda x: x):
    return bisect.bisect_left(sorted_list, target, key=key)

def is_prevalent(candidate, instance_counts, bounding_boxes, minprev):
    feature_counts = {feature: instance_counts[feature] for feature in candidate}
    if len(feature_counts) != len(candidate):
        return False

    participation_counts = defaultdict(int)
    for feature in candidate:
        sorted_instances = bounding_boxes[feature]
        for other_feature in candidate - {feature}:
            count = 0
            for bb in bounding_boxes[other_feature]:
                index = binary_search(sorted_instances, bb.aabb_min[0], key=lambda x: x.aabb_max[0])
                count += sum(1 for inst in sorted_instances[index:] if inst.intersects(bb))
            participation_counts[feature] = max(participation_counts[feature], count)

    prevalences = [participation_counts[feature] / feature_counts[feature] for feature in candidate]
    return min(prevalences) >= minprev if prevalences else False

def process_time_step(time, state, maxdist, distance_metric, minprev, angle, all_features, prevalent_pairs, time_prevalent_pairs, bounding_boxes_per_time):
    current_features = set(inst.id.feature for inst in state.instances)
    all_features.update(current_features)
    instance_counts = state.count_instances()
    
    bounding_boxes = defaultdict(list)
    for inst in state.instances:
        bb = create_rotated_bounding_box(inst, maxdist, distance_metric, angle)
        bounding_boxes[inst.id.feature].append(bb)
    
    for feature in bounding_boxes:
        bounding_boxes[feature].sort(key=lambda bb: bb.aabb_min[0])
    
    bounding_boxes_per_time[time] = bounding_boxes

    sorted_features = sorted(current_features)
    for i, f1 in enumerate(sorted_features):
        for f2 in sorted_features[i+1:]:
            pair = frozenset([f1, f2])
            if is_prevalent(pair, instance_counts, bounding_boxes, minprev):
                prevalent_pairs[time].add(pair)
                time_prevalent_pairs.add(pair)

def step1_build_bounding_boxes_and_find_prevalent_pairs(input_data, maxdist, distance_metric, minprev, angle=45):
    all_features = set()
    bounding_boxes_per_time = {}
    prevalent_pairs = defaultdict(set)
    time_prevalent_pairs = set()

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_time_step, time, state, maxdist, distance_metric, minprev, angle, all_features, prevalent_pairs, time_prevalent_pairs, bounding_boxes_per_time) for time, state in input_data.states.items()]
        for future in concurrent.futures.as_completed(futures):
            future.result()

    return all_features, bounding_boxes_per_time, prevalent_pairs, time_prevalent_pairs

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

def step3_and_4_generate_and_prune_candidates(input_data, time_prevalent_pairs, bounding_boxes_per_time, minfreq, minprev, prevalent_pairs):
    candidates = list(time_prevalent_pairs)
    prevalence_cache = {}
    result = list(time_prevalent_pairs)  # Start with prevalent pairs
    k = 3  # Start from size 3 candidates

    while candidates:
        print(f"Generating candidates of size {k}")
        new_candidates = generate_candidates(candidates, k)
        prevalent_candidates = []

        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_to_candidate = {executor.submit(check_prevalence, candidate, bounding_boxes_per_time, input_data, minfreq, minprev, prevalence_cache, prevalent_pairs): candidate for candidate in new_candidates}
            for future in concurrent.futures.as_completed(future_to_candidate):
                candidate = future_to_candidate[future]
                is_prevalent = future.result()
                if is_prevalent:
                    prevalent_candidates.append(candidate)

        result.extend(prevalent_candidates)
        candidates = prevalent_candidates
        print(f"Found {len(prevalent_candidates)} prevalent candidates of size {k}")
        k += 1
    
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
            print(f"Feature {feature} not found in instance_counts.")
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

    # Step 1
    step1_start = tm.time()
    all_features, bounding_boxes_per_time, prevalent_pairs, time_prevalent_pairs = step1_build_bounding_boxes_and_find_prevalent_pairs(input_data, maxdist, distance_metric, minprev)
    step1_end = tm.time()
    if verbose > 0:
        print(f"Step 1 (Build bounding boxes and find prevalent pairs) took {step1_end - step1_start:.2f} seconds")

    # Step 2
    step2_start = tm.time()
    time_prevalent_pairs = step2_find_time_prevalent_pairs(prevalent_pairs, time_prevalent_pairs, minfreq, len(input_data.states))
    step2_end = tm.time()
    if verbose > 0:
        print(f"Step 2 (Find time-prevalent pairs) took {step2_end - step2_start:.2f} seconds")

    if not time_prevalent_pairs:
        return []

    # Step 3 & 4
    step34_start = tm.time()
    result = step3_and_4_generate_and_prune_candidates(input_data, time_prevalent_pairs, bounding_boxes_per_time, minfreq, minprev, prevalent_pairs)
    step34_end = tm.time()
    if verbose > 0:
        print(f"Steps 3 & 4 (Generate and prune candidates) took {step34_end - step34_start:.2f} seconds")

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
