import time as tm
from collections import defaultdict
from itertools import combinations
from algorithms.utils import manhattan_distance, chebyshev_distance

class BoundingBox:
    def __init__(self, min_coords, max_coords):
        self.min_coords = min_coords
        self.max_coords = max_coords

    def intersects(self, other):
        return all(self.min_coords[i] <= other.max_coords[i] and self.max_coords[i] >= other.min_coords[i] 
                   for i in range(len(self.min_coords)))

def create_bounding_box(instance, distance, distance_func):
    if distance_func == manhattan_distance:
        offset = distance / 2
    elif distance_func == chebyshev_distance:
        offset = distance
    else:
        raise ValueError("Unsupported distance function")
    
    min_coords = [instance.pos.x - offset, instance.pos.y - offset]
    max_coords = [instance.pos.x + offset, instance.pos.y + offset]
    return BoundingBox(min_coords, max_coords)

def is_prevalent(candidate, instance_counts, bounding_boxes, minprev):
    feature_counts = {}
    for feature in candidate:
        if feature not in instance_counts or feature not in bounding_boxes:
            return False
        feature_counts[feature] = instance_counts[feature]

    participation_counts = defaultdict(int)
    for feature in candidate:
        for other_feature in candidate - {feature}:
            count = sum(1 for bb1 in bounding_boxes[feature] 
                        for bb2 in bounding_boxes[other_feature] 
                        if bb1.intersects(bb2))
            participation_counts[feature] = max(participation_counts[feature], count)

    prevalences = [participation_counts[feature] / feature_counts[feature] for feature in candidate]
    return min(prevalences) >= minprev if prevalences else False

def generate_candidates(prev_candidates, k):
    new_candidates = set()
    for i, c1 in enumerate(prev_candidates):
        for c2 in prev_candidates[i+1:]:
            new_candidate = c1.union(c2)
            if len(new_candidate) == k:
                if all(frozenset(comb) in prev_candidates for comb in combinations(new_candidate, k-1)):
                    new_candidates.add(new_candidate)
    return new_candidates

def BBmaxspatiotempcolloc(input_data, maxdist, minprev, minfreq, distance_metric=chebyshev_distance, predictions=False, verbose=0, clean_trees=False):
    start_time = tm.time()
    result = []
    all_features = set()
    bounding_boxes_per_time = {}

    # Step 1: Build bounding boxes and find prevalent pairs for each time point
    prevalent_pairs = defaultdict(set)
    step1_start = tm.time()

    for time, state in input_data.states.items():
        current_features = set(inst.id.feature for inst in state.instances)
        all_features.update(current_features)
        instance_counts = state.count_instances()
        
        bounding_boxes = defaultdict(list)
        for inst in state.instances:
            bb = create_bounding_box(inst, maxdist, distance_metric)
            bounding_boxes[inst.id.feature].append(bb)
        bounding_boxes_per_time[time] = bounding_boxes
        
        for f1 in current_features:
            for f2 in current_features:
                if f1 >= f2:
                    continue
                pair = frozenset([f1, f2])
                if is_prevalent(pair, instance_counts, bounding_boxes, minprev):
                    prevalent_pairs[time].add(pair)

    step1_end = tm.time()
    if verbose > 0:
        print(f"Step 1 (Build bounding boxes and find prevalent pairs) took {step1_end - step1_start:.2f} seconds")

    # Step 2: Find time-prevalent pairs
    step2_start = tm.time()
    time_prevalent_pairs = {
        pair for pair in set.union(*prevalent_pairs.values())
        if sum(pair in pairs for pairs in prevalent_pairs.values()) >= minfreq * len(input_data.states)
    }
    step2_end = tm.time()
    if verbose > 0:
        print(f"Step 2 (Find time-prevalent pairs) took {step2_end - step2_start:.2f} seconds")

    if not time_prevalent_pairs:
        return []

    # Step 3 & 4: Generate and prune candidate set
    candidates = list(time_prevalent_pairs)
    k = 3

    step34_start = tm.time()
    while candidates:
        new_candidates = generate_candidates(candidates, k)
        prevalent_candidates = []

        for candidate in new_candidates:
            time_prevalent_count = 0
            for time, state in input_data.states.items():
                if all(feature in state.count_instances() for feature in candidate):
                    instance_counts = state.count_instances()
                    if is_prevalent(candidate, instance_counts, bounding_boxes_per_time[time], minprev):
                        time_prevalent_count += 1
                        if time_prevalent_count >= minfreq * len(input_data.states):
                            prevalent_candidates.append(candidate)
                            break

        result.extend(prevalent_candidates)
        candidates = prevalent_candidates
        k += 1

        if verbose > 1:
            print(f"Iteration for k={k-1}: {len(prevalent_candidates)} prevalent candidates found")

    step34_end = tm.time()
    if verbose > 0:
        print(f"Steps 3 & 4 (Generate and prune candidates) took {step34_end - step34_start:.2f} seconds")

    # Ensure maximality
    step5_start = tm.time()
    maximal_result = []
    for pattern in sorted(result, key=len, reverse=True):
        if not any(pattern.issubset(maximal) for maximal in maximal_result):
            maximal_result.append(pattern)
    step5_end = tm.time()
    if verbose > 0:
        print(f"Step 5 (Ensure maximality) took {step5_end - step5_start:.2f} seconds")

    end_time = tm.time()
    if verbose > 0:
        print(f"Total execution time: {end_time - start_time:.2f} seconds")
    
    return maximal_result