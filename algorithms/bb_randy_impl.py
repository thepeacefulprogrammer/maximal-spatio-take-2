import itertools
from time import perf_counter
import numpy as np
from algorithms.planesweep import planesweep
from structures.graph import Graph
import bisect

def manhattan_distance(p1, p2):
    return sum(abs(a - b) for a, b in zip(p1, p2))

def chebyshev_distance(p1, p2):
    return max(abs(a - b) for a, b in zip(p1, p2))

def bounding_box(points):
    min_point = tuple(map(min, zip(*points)))
    max_point = tuple(map(max, zip(*points)))
    return min_point, max_point

def object_bounding_box(position, d, distance_metric):
    q1, q2 = position.x, position.y
    min_point = (q1 - d/2, q2 - d/2)
    max_point = (q1 + d/2, q2 + d/2)
    return min_point, max_point

def boxes_intersect(bb1, bb2):
    (min1, max1), (min2, max2) = bb1, bb2
    return (min1[0] <= max2[0] and max1[0] >= min2[0] and
            min1[1] <= max2[1] and max1[1] >= min2[1])

def composite_bounding_box(bounding_boxes, maxdist, distance_metric):
    if not bounding_boxes:
        return None
    min_points, max_points = zip(*bounding_boxes)
    return (
        (min(p[0] for p in min_points), min(p[1] for p in min_points)),
        (max(p[0] for p in max_points), max(p[1] for p in max_points))
    )

def find_prevalent_pairs(input, maxdist, minprev):
    all_pairs = set()
    for state in input.states.values():
        state.build_local_ICPI(maxdist)
        state.find_prevalent_pairs(minprev)
        all_pairs.update(state.prevalent_pairs)
    return all_pairs

def find_nonfrequent_pairs(input, all_pairs, minfreq_it):
    forbidden_pairs = set()
    for pair in all_pairs:
        freq_it = sum(1 for x in input.states.values() if pair in x.prevalent_pairs)
        if freq_it < minfreq_it:
            forbidden_pairs.add(frozenset(pair))
    return forbidden_pairs

def remove_nonfrequent_pairs(input, forbidden_pairs):
    for state in input.states.values():
        state.prevalent_pairs = [pair for pair in state.prevalent_pairs if frozenset(pair) not in forbidden_pairs]

def build_local_and_global_candidate_sets(input):
    initial_global_candidates = {}
    for state in input.states.values():
        state.build_local_candidates()
        if state.local_candidates:
            for c in state.local_candidates:
                if len(c) not in initial_global_candidates:
                    initial_global_candidates[len(c)] = set()
                initial_global_candidates[len(c)].add(frozenset(c))
    return initial_global_candidates

def prune_candidate_set(input, initial_global_candidates, minfreq_it):
    global_candidates = {}
    while len(initial_global_candidates) > 0:
        maxlen = max(initial_global_candidates.keys())
        global_candidates[maxlen] = set()
        candidates = initial_global_candidates.pop(maxlen)
        
        for candidate in candidates:
            freq_it = sum(1 for x in input.states.values() if any(candidate.issubset(a) for a in x.local_candidates))
            
            if freq_it >= minfreq_it:
                global_candidates[maxlen].add(candidate)
            else:
                if len(candidate) > 2:
                    for item in candidate:
                        if maxlen - 1 not in initial_global_candidates:
                            initial_global_candidates[maxlen - 1] = set()
                        initial_global_candidates[maxlen - 1].add(frozenset(candidate.difference({item})))
        
        if len(global_candidates[maxlen]) == 0:
            global_candidates.pop(maxlen)
        
        if maxlen == 2:
            break
    
    return global_candidates

def BBmaxspatiotempcolloc(input, maxdist, minprev, minfreq, distance_metric=manhattan_distance, predictions=False, verbose=0, clean_trees=0):
    result = []
    minfreq_it = minfreq * len(input.states)

    all_pairs = find_prevalent_pairs(input, maxdist, minprev)
    if not all_pairs:
        print("No prevalent pairs found.")
        return result

    forbidden_pairs = find_nonfrequent_pairs(input, all_pairs, minfreq_it)
    remove_nonfrequent_pairs(input, forbidden_pairs)

    initial_global_candidates = build_local_and_global_candidate_sets(input)
    if not initial_global_candidates:
        print("No global candidates found.")
        return result

    global_candidates = prune_candidate_set(input, initial_global_candidates, minfreq_it)

    while len(global_candidates) > 0:
        maxlen = max(global_candidates.keys())

        if maxlen == 2:
            result.extend(all_pairs)
            break
        else:
            candidates = global_candidates.pop(maxlen)
            new_candidates = set()

            for candidate in candidates:
                frequency = []
                for statenum, state in input.states.items():
                    candidate_features = list(candidate)
                    
                    # Get all instances for each feature in the candidate
                    candidate_instances = [
                        instance for feature in candidate_features
                        for instance in state.instances if instance.id.feature == feature
                    ]
                    
                    if not candidate_instances:
                        continue  # Skip this state if no instances found
                    
                    bounding_boxes = [object_bounding_box(instance.pos, maxdist, distance_metric) for instance in candidate_instances]
                    Bc = composite_bounding_box(bounding_boxes, maxdist, distance_metric)
                    
                    if Bc is None:
                        continue  # Skip if we couldn't create a bounding box
                    
                    Ic = candidate_instances
                    for j in range(2):  # Assuming 2D space
                        Ic, Mc = improved_instance_identification(state, Bc, maxdist, distance_metric, j)
                        if not Ic:
                            break
                        Bc = Mc
                    
                    if Ic:
                        # Calculate prevalence based on the number of instances for each feature
                        prevalences = []
                        for feature in candidate_features:
                            feature_instances = [i for i in Ic if i.id.feature == feature]
                            total_feature_instances = len([i for i in state.instances if i.id.feature == feature])
                            if total_feature_instances > 0:
                                prevalences.append(len(feature_instances) / total_feature_instances)
                        
                        if prevalences and min(prevalences) >= minprev:
                            frequency.append(statenum)

                if len(frequency) >= minfreq_it:
                    result.append(candidate)
                else:
                    for item in candidate:
                        new_candidates.add(frozenset(candidate.difference({item})))

            if len(new_candidates) > 0:
                if maxlen - 1 not in global_candidates:
                    global_candidates[maxlen - 1] = set()
                global_candidates[maxlen - 1].update(new_candidates)

    return result

def improved_instance_identification(state, Bc, maxdist, distance_metric, j):
    Ic = []
    Mc = Bc

    T = sorted(state.instances, key=lambda x: object_bounding_box(x.pos, maxdist, distance_metric)[0][j])
    k = bisect.bisect_left(T, Bc[0][j], key=lambda x: object_bounding_box(x.pos, maxdist, distance_metric)[1][j])

    while k < len(T) and object_bounding_box(T[k].pos, maxdist, distance_metric)[0][j] <= Bc[1][j]:
        instance_bb = object_bounding_box(T[k].pos, maxdist, distance_metric)
        if boxes_intersect(Mc, instance_bb):
            Ic.append(T[k])
            Mc = composite_bounding_box([Mc, instance_bb], maxdist, distance_metric)
        k += 1

    return Ic, Mc