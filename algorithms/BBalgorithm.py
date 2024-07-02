from functools import reduce
from time import perf_counter
import numpy as np

# Calculate Manhattan distance (dist1)
def manhattan_distance(p1, p2):
    return sum(abs(a - b) for a, b in zip(p1, p2))

# Calculate Chebyshev distance (dist∞)
def chebyshev_distance(p1, p2):
    return max(abs(a - b) for a, b in zip(p1, p2))

# Calculate bounding box for a set of points
def bounding_box(points):
    min_point = tuple(map(min, zip(*points)))
    max_point = tuple(map(max, zip(*points)))
    return min_point, max_point

# Calculate the bounding box for an object with given position and distance threshold using Manhattan distance
def object_bounding_box_manhattan(position, d):
    q1, q2 = position
    min_point = (q1 - d, q2 - d)
    max_point = (q1 + d, q2 + d)
    return min_point, max_point

# Calculate the bounding box for an object with given position and distance threshold using Chebyshev distance
def object_bounding_box_chebyshev(position, d):
    q1, q2 = position
    min_point = (q1 - d, q2 - d)
    max_point = (q1 + d, q2 + d)
    return min_point, max_point

# Rotate points by π/4 radians
def rotate_point(point):
    q1, q2 = point
    rotated_q1 = (q1 - q2) / np.sqrt(2)
    rotated_q2 = (q1 + q2) / np.sqrt(2)
    return rotated_q1, rotated_q2

# Rotate bounding box
def rotate_bounding_box(min_point, max_point):
    rotated_min = rotate_point(min_point)
    rotated_max = rotate_point(max_point)
    return bounding_box([rotated_min, rotated_max])

# Check if two bounding boxes intersect
def boxes_intersect(bb1, bb2):
    (min1, max1), (min2, max2) = bb1, bb2
    return (min1[0] <= max2[0] and max1[0] >= min2[0] and
            min1[1] <= max2[1] and max1[1] >= min2[1])

# Create Composite Bounding Box (CBB) for a set of points
def composite_bounding_box(points, d, distance_fn):
    rotated_points = [rotate_point(p) for p in points]
    bounding_boxes = [object_bounding_box_manhattan(p, d) if distance_fn == manhattan_distance else object_bounding_box_chebyshev(p, d) for p in rotated_points]
    min_points, max_points = zip(*bounding_boxes)
    return bounding_box(min_points), bounding_box(max_points)
def output_prevalent_pairs(input, all_pairs, verbose):
    if verbose > 1:
        print("Prevalent pairs for each time moment:")
        for time in input.states.keys():
            print("Time t=", time)
            state = input.states[time]
            print(state.prevalent_pairs)
        print("All pairs:")
        print(all_pairs)


def output_nonfrequent_pairs(forbidden_pairs, verbose):
    if verbose > 1:
        print("Non frequent pairs")
        print(forbidden_pairs)


def output_size2_MDCOPS(input, verbose):
    if verbose > 1:
        print("Prevalent and frequent pairs for each time moment:")
        for time in input.states.keys():
            print("Time t=", time)
            state = input.states[time]
            print(state.prevalent_pairs)


def output_global_candidates(global_candidates, verbose):
    if verbose > 1:
        print("Pruned global candidates:")
        print(global_candidates)


def output_local_and_global_candidate_sets(input, initial_global_candidates, verbose):
    if verbose > 1:
        print("Local candidates:")
        for time in input.states.keys():
            print("Time t=", time)
            state = input.states[time]
            print(state.local_candidates)

        print("Initial global candidates (including non frequent ones)")
        print(initial_global_candidates)


def output_results_length_cache(maxlen, result, cache, verbose):
    if verbose > 1:
        print("Current length ", maxlen)
        print("\tCurrent result set: ", result)
        print("\tCurrent cache: ", cache)


def output_step_time(step, time, verbose):
    if verbose > 0:
        print("Step", step, "time", time)


def output_iteration_time(maxlen, time, verbose):
    if verbose > 0:
        print("Main, length ", maxlen, ":", time)


def output_cache_stats(cache, cachetime, cachehit, cachemiss, verbose):
    if verbose > 0:
        print("Cache time", cachetime, "cache hit ", cachehit, "cachemiss ", cachemiss, "cachesize",
              reduce(lambda x, y: x + y, (len(x) for x in cache.values())))


def output_glocal_candidate_stats(maxlen, global_candidates, verbose):
    if verbose > 0:
        if maxlen - 1 in global_candidates:
            print("Candidate size ", maxlen - 1, " count ", len(global_candidates[maxlen - 1]))
        else:
            print("Candidate size ", maxlen - 1, " count 0")


def output_iteration_state(maxlen, result, cache, new_candidates, verbose):
    if verbose > 1:
        print("Current length ", maxlen)
        print("\tCurrent result set: ", result)
        print("\tCurrent cache: ", cache)
        print("\tNew candidates added:", new_candidates)


# Finds prevalent pairs in each time moment (stored as a field prevalent_pairs in the WorldState objects)
# as well as a set of all prevalent pairs (returned as a function value)
def find_prevalent_pairs(input, maxdist, minprev):
    # Declare a set which will contain all frequent pairs
    all_pairs = set()

    # Find prevalent pairs in each time moment and put them all into one set
    for state in input.states.values():
        # Build an iCPI tree for the time moment represented by a state
        # iCPI tree is built using plane sweep
        state.build_local_ICPI(maxdist)
        # Convert neighborhood information from the iCPI tree into a set of prevalent pairs (local variable in the state object) 
        state.find_prevalent_pairs(minprev)
        # Compute set sum of the up-to-date found pairs with the pairs found for the current state
        # Note that duplicates are removed since all_pairs is a set
        all_pairs.update(state.prevalent_pairs)

    # Return the found set of frequent pairs
    return all_pairs


# Finds pairs that are prevalent in some time moments, but are infrequent (returned as a function value)
def find_nonfrequent_pairs(input, all_pairs, minfreq_it):
    # Declare a set for storing infrequent pairs
    forbidden_pairs = set()

    # For all pair prevalent at some time moments 
    for pair in all_pairs:
        # Compute frequency of a pair expressed as a number of time moments in which it is prevalent
        freq_it = sum(1 for x in input.states.values() if pair in x.prevalent_pairs)
        # If the computed frequency is smaller than the minfreq threshold then add the pair to the result set
        if freq_it < minfreq_it:
            forbidden_pairs.add(pair)

    # Return the found set of infrequent pairs
    return forbidden_pairs


# Modifies states of each time moment by removing nonfrequent pairs from prevalent_pairs sets
def remove_nonfrequent_pairs(input, forbidden_pairs):
    # Update prevalent_pairs in all time moments to store only time prevalent pairs
    for state in input.states.values():
        state.remove_prevalent_pairs(forbidden_pairs)


# Generates local candidate sets and their union: a global candidate set
def build_local_and_global_candidate_sets(input):
    # Declare a map for storing global candidates
    # Each key will represent size of a candidate set. The corresponding value will be a set of candidate sets of proper size 
    initial_global_candidates = {}

    # For each time moment...
    for state in input.states.values():
        # ...build a set of local candidates based on frequent prevalent pairs
        state.build_local_candidates()
        # For each of the local candidates insert it into an appropriate set in the initial_global_candidates map 
        for c in state.local_candidates:
            # If the appropriate entry in the initial_global_candidates map is not created yet then create it
            if not len(c) in initial_global_candidates:
                initial_global_candidates[len(c)] = set()
            # Insert the candidate c into appropriate set of candidates
            initial_global_candidates[len(c)].update({c})

    # Return the set of global_candidates
    return initial_global_candidates


# Modifies initial global candidate set by removing candidates that are infrequent.
# This is done by checking whether a candidate is a subset of enough local candidates to be frequent
# In case a candidate is not frequent among other candidates, it is replaced by all of its size-1 subsets   
def prune_candidate_set(input, initial_global_candidates, minfreq_it):
    # Declare map for storing pruned candidate sets
    # Each key will represent size of a candidate set. The corresponding value will be a set of candidate sets of proper size 
    global_candidates = {}

    # As long as there are still sets of candidates in the initial global candidates map...
    while len(initial_global_candidates) > 0:
        # Get maximal candidate length yet unprocessed
        maxlen = max(initial_global_candidates.keys())
        # Initialize an entry in the global_candidates map for candidates of size maxlen
        global_candidates[maxlen] = set()

        # Extract and remove largest candidates from the initial global candidates map
        candidates = initial_global_candidates[maxlen]
        initial_global_candidates.pop(maxlen)

        # For each of the extracted candidates
        for candidate in candidates:
            # Check if superset is already added, if so, omit the candidate
            remove = False
            for k in global_candidates.keys():
                if any(candidate <= a for a in global_candidates[k] if k > maxlen):
                    remove = True
                    break
            if remove:
                continue

            # Compute frequency of a candidate
            freq_it = reduce(lambda x, y: x + y,
                             (1 for x in input.states.values() if any(candidate <= a for a in x.local_candidates)))

            # If candidate is potentially frequent then add it to filtered global candidates
            # Otherwise test all of its subsets by adding them to the initial global candidates         
            if freq_it >= minfreq_it:
                global_candidates[maxlen].add(candidate)
            else:
                if len(candidate) > 2:  # Do not add subsets of size 2 sets
                    for item in candidate:
                        if not maxlen - 1 in initial_global_candidates:
                            initial_global_candidates[maxlen - 1] = set()
                        initial_global_candidates[maxlen - 1].add(candidate.difference({item}))

        # If no global candidates of length maxlen are potentially frequent then remove the corresponding entry from the global candidates map
        if len(global_candidates[maxlen]) == 0:
            global_candidates.pop(maxlen)

        # If maximal candidate length equals 2 then break the loop
        # Notice that in the first step we are always processing size-2 candidates. Hence the algorithm never terminates at this point
        if maxlen == 2:
            break

    # Return the pruned set of global candidates
    return global_candidates


def find_frequency_without_estimations(input, candidate, cache, minfreq_it, minprev, frequency, from_cache, save_trees):
    cachetime = 0;
    cachehit = 0;
    cachemiss = 0;

    for statenum in input.states:  # for each of the time moments
        state = input.states[statenum]

        # check in cache if the candidate is a subset of a locally  prevalent candidate
        t3 = perf_counter()
        incache = any(candidate <= a for a in cache[statenum])
        t4 = perf_counter()
        cachetime = cachetime + t4 - t3

        if incache:
            cachehit = cachehit + 1
            frequency.append(statenum)  # remember the time moment at which the candidate is prevalent
            from_cache.add(
                statenum)  # remember the time moment at which the cache allowed to determine that the candidate is prevalent
        else:
            cachemiss = cachemiss + 1
            # if not cached, then compute candidate's prevalence in the current time moment
            prevalence = state.compute_prevalence(candidate, save_trees);
            if prevalence >= minprev:
                frequency.append(statenum)  # remember the time moment at which the candidate is  prevalent

    return (cachetime, cachehit, cachemiss)


def clean_tree_memory(input,maxlen):
    for statenum in input.states: # for each of the time moments
        state=input.states[statenum]
        state.clean_memory(maxlen) #delete every tree level at maxlen depth


def update_global_candidates(global_candidates, new_candidates, maxlen, result):
    # check if new candidates are not subsets of current results
    for candidate in new_candidates:
        if not any(candidate <= a for a in result):  # If new candidate is not a subset of any result then...
            if not maxlen - 1 in global_candidates:
                global_candidates[maxlen - 1] = set()
            global_candidates[maxlen - 1].add(candidate)  # ... add it to global_candidates


# Finds maximal spatial-temporal collocations
def BBmaxspatiotempcolloc(input, maxdist, minprev, minfreq, distance_metric=manhattan_distance, predictions=False,
                        save_trees=0, verbose=0, clean_trees=0):
    # Prepare result list and compute minimal number of time moments to qualify a pair as frequent
    result = []
    minfreq_it = minfreq * len(input.states)

    # Measure the initial time to compute prevalent pairs
    start = perf_counter()
    all_pairs = find_prevalent_pairs(input, maxdist, minprev)
    output_prevalent_pairs(input, all_pairs, verbose)

    # Measure time for the first step
    end = perf_counter()
    output_step_time(1, end - start, verbose)

    # Measure the initial time to compute non-frequent pairs
    start = perf_counter()
    forbidden_pairs = find_nonfrequent_pairs(input, all_pairs, minfreq_it)
    output_nonfrequent_pairs(forbidden_pairs, verbose)

    # Remove the non-frequent pairs
    remove_nonfrequent_pairs(input, forbidden_pairs)
    output_size2_MDCOPS(input, verbose)

    # Measure time for the second step
    end = perf_counter()
    output_step_time(2, end - start, verbose)

    # Measure the initial time to compute local and global candidate sets
    start = perf_counter()
    initial_global_candidates = build_local_and_global_candidate_sets(input)
    output_local_and_global_candidate_sets(input, initial_global_candidates, verbose)

    # Measure time for the third step
    end = perf_counter()
    output_step_time(3, end - start, verbose)

    # Measure the initial time to prune the candidate set
    start = perf_counter()
    global_candidates = prune_candidate_set(input, initial_global_candidates, minfreq_it)
    output_global_candidates(global_candidates, verbose)

    # Measure time for the fourth step
    end = perf_counter()
    output_step_time(4, end - start, verbose)

    # Initiaze cache map storing already computed results
    cache = {x: set() for x in input.states.keys()}

    # While there are still candidates to be processed...
    while len(global_candidates) > 0:
        # Measure the initial time to process the next iteration
        start = perf_counter()

        # Get maximal candidate length yet unprocessed
        maxlen = max(global_candidates.keys())

        # If clean_trees and save_trees options are enabled then clean tree memory
        if clean_trees and save_trees:
            clean_tree_memory(input, maxlen)

        # If maximal candidate length equals 2 then process size-2 candidates
        if maxlen == 2:
            # All previously computed frequent pairs are stored in all_pairs variable
            # Forbidden pairs are those that were non-frequent
            # All of the frequent and prevalent pairs are stored in the result list
            result.extend(all_pairs)

            # Output current state of the algorithm
            output_results_length_cache(maxlen, result, cache, verbose)
        else:
            # Otherwise process candidates larger than size-2
            # Extract and remove largest candidates from the global candidates map
            candidates = global_candidates.pop(maxlen)

            # Initialize set for storing new candidate sets
            new_candidates = set()
            # Initialize cache statistics counters
            cachetime = 0
            cachehit = 0
            cachemiss = 0

            # For each of the candidates...
            for candidate in candidates:
                # Initialize variables for storing frequencies, candidates taken from cache and total time taken by cache operations
                frequency = []
                from_cache = set()

                # Compute candidate frequencies in a list of time moments in which the candidate is frequent
                # Also update cache with found candidate sets
                # Returns cache statistics tuple (cachetime, cachehit, cachemiss)
                cache_stats = find_frequency_without_estimations(input, candidate, cache, minfreq_it, minprev,
                                                                 frequency, from_cache, save_trees)
                cachetime += cache_stats[0]
                cachehit += cache_stats[1]
                cachemiss += cache_stats[2]

                # If candidate is frequent then add it to the result set
                if len(frequency) >= minfreq_it:
                    result.append(candidate)
                # Otherwise add it to the set of new candidates
                else:
                    for item in candidate:
                        new_candidates.add(candidate.difference({item}))

            # If there are new candidates, then update global candidate set
            if len(new_candidates) > 0:
                update_global_candidates(global_candidates, new_candidates, maxlen, result)

            # Output cache statistics
            output_cache_stats(cache, cachetime, cachehit, cachemiss, verbose)

            # Output candidate statistics
            output_glocal_candidate_stats(maxlen, global_candidates, verbose)

            # Output current state of the algorithm
            output_iteration_state(maxlen, result, cache, new_candidates, verbose)

        # Measure time for the current iteration
        end = perf_counter()
        output_iteration_time(maxlen, end - start, verbose)

    # Return the final result
    return result

# Usage example
# Make sure to define the WorldState class and other helper functions as needed, and instantiate 'input' appropriately
