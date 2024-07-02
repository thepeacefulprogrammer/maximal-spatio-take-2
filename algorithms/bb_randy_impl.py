from functools import reduce
from time import perf_counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
from algorithms.mainalg import (
    output_prevalent_pairs,
    output_nonfrequent_pairs,
    output_size2_MDCOPS,
    output_global_candidates,
    output_local_and_global_candidate_sets,
    output_results_length_cache,
    output_step_time,
    output_iteration_time,
    output_cache_stats,
    output_glocal_candidate_stats,
    output_iteration_state,
    find_prevalent_pairs,
    find_nonfrequent_pairs,
    remove_nonfrequent_pairs,
    build_local_and_global_candidate_sets,
    prune_candidate_set,
    find_frequency_without_estimations,
    clean_tree_memory,
    update_global_candidates
)

from structures.input import WorldState

# Calculate Manhattan distance (dist1)
def manhattan_distance(p1, p2):
    return sum(abs(a - b) for a, b in zip(p1, p2))

def binary_search_first_neighbor(T_c, o1, intL, M_ic_prime):
    neighbors = T_c[o1]
    low, high = 0, len(neighbors) - 1

    while low <= high:
        mid = (low + high) // 2
        if intL(M_ic_prime) <= neighbors[mid].intU:
            high = mid - 1
        else:
            low = mid + 1

    return low
    
def intU(M):
    # Implementation for intU (upper bound of the interval)
    return max(M[1])

def intL(M):
    # Implementation for intL (lower bound of the interval)
    return min(M[0])

def intersection(bb1, bb2):
    (min1, max1), (min2, max2) = bb1, bb2
    min_point = (max(min1[0], min2[0]), max(min1[1], min2[1]))
    max_point = (min(max1[0], max2[0]), min(max1[1], max2[1]))

    if min_point[0] <= max_point[0] and min_point[1] <= max_point[1]:
        return (min_point, max_point)
    else:
        return None
        
# Improved instance identification with binary search
def improved_instance_identification(candidate, prefix_instances, prefix_CBBs, icpi_tree, sorting_dim):
    candidate_instances = []
    candidate_bounding_boxes = []

    for k, (ic_prime, M_ic_prime) in enumerate(zip(prefix_instances, prefix_CBBs)):
        if len(ic_prime) == 0:
            continue
        
        o1 = next(iter(ic_prime))  # Use next(iter()) for clarity and efficiency
        
        k = binary_search_first_neighbor(icpi_tree.T_c, o1, intL, M_ic_prime)
        
        while k < len(icpi_tree.T_c[o1]):
            oc = icpi_tree.T_c[o1][k]
            if intU(M_ic_prime) > intL(icpi_tree.CBB[oc]):
                break
            
            M_ic = intersection(M_ic_prime, icpi_tree.CBB[oc])
            if M_ic is not None:
                ic = ic_prime | {oc}
                candidate_instances.append(ic)
                candidate_bounding_boxes.append(M_ic)
            k += 1
    
    return candidate_instances, candidate_bounding_boxes

def process_candidates_parallel(candidates, prefix_instances, prefix_CBBs, icpi_tree, sorting_dim):
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(improved_instance_identification, candidate, prefix_instances, prefix_CBBs, icpi_tree, sorting_dim)
                   for candidate in candidates]
        results = [future.result() for future in as_completed(futures)]
    return results

def BBmaxspatiotempcolloc(input, maxdist, minprev, minfreq, distance_metric=manhattan_distance, predictions=False, save_trees=0, verbose=0, clean_trees=0):
    # Prepare result list and compute minimal number of time moments to qualify a pair as frequent
    result = []
    minfreq_it = minfreq * len(input.states)

    # Measure the initial time to compute prevalent pairs
    start = perf_counter()
    all_pairs = find_prevalent_pairs(input, maxdist, minprev)
    output_prevalent_pairs(input, all_pairs, verbose)
    end = perf_counter()
    output_step_time(1, end - start, verbose)

    # Measure the initial time to compute non-frequent pairs
    start = perf_counter()
    forbidden_pairs = find_nonfrequent_pairs(input, all_pairs, minfreq_it)
    output_nonfrequent_pairs(forbidden_pairs, verbose)
    remove_nonfrequent_pairs(input, forbidden_pairs)
    output_size2_MDCOPS(input, verbose)
    end = perf_counter()
    output_step_time(2, end - start, verbose)

    # Measure the initial time to compute local and global candidate sets
    start = perf_counter()
    initial_global_candidates = build_local_and_global_candidate_sets(input)
    output_local_and_global_candidate_sets(input, initial_global_candidates, verbose)
    end = perf_counter()
    output_step_time(3, end - start, verbose)

    # Measure the initial time to prune the candidate set
    start = perf_counter()
    global_candidates = prune_candidate_set(input, initial_global_candidates, minfreq_it)
    output_global_candidates(global_candidates, verbose)
    end = perf_counter()
    output_step_time(4, end - start, verbose)

    # Initialize cache map storing already computed results
    cache = {x: set() for x in input.states.keys()}

    # Main loop to process candidates
    while len(global_candidates) > 0:
        start = perf_counter()

        # Get maximal candidate length yet unprocessed
        maxlen = max(global_candidates.keys())

        if clean_trees and save_trees:
            clean_tree_memory(input, maxlen)

        if maxlen == 2:
            result.extend(all_pairs)
            output_results_length_cache(maxlen, result, cache, verbose)
        else:
            candidates = global_candidates.pop(maxlen)
            new_candidates = set()
            cachetime, cachehit, cachemiss = 0, 0, 0

            for candidate in candidates:
                prefix_instances = [[] for _ in input.states]
                prefix_CBBs = [[] for _ in input.states]
                icpi_tree = list(input.states.values())[0].icpi
                sorting_dim = 0

                candidate_instances, candidate_bounding_boxes = process_candidates_parallel(
                    [candidate], prefix_instances, prefix_CBBs, icpi_tree, sorting_dim
                )[0]  # Unpack the single result tuple

                frequency = []
                from_cache = set()
                cache_stats = find_frequency_without_estimations(
                    input, candidate, cache, minfreq_it, minprev, frequency, from_cache, save_trees
                )
                cachetime += cache_stats[0]
                cachehit += cache_stats[1]
                cachemiss += cache_stats[2]

                if len(frequency) >= minfreq_it:
                    result.append(candidate)
                else:
                    new_candidates.update(candidate.difference({item}) for item in candidate)

            if new_candidates:
                update_global_candidates(global_candidates, new_candidates, maxlen, result)

            output_cache_stats(cache, cachetime, cachehit, cachemiss, verbose)
            output_glocal_candidate_stats(maxlen, global_candidates, verbose)
            output_iteration_state(maxlen, result, cache, new_candidates, verbose)

        end = perf_counter()
        output_iteration_time(maxlen, end - start, verbose)

    return result
