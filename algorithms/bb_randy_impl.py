import math
from collections import defaultdict
from itertools import combinations

class Position:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

class FeatureInstanceIdentifier:
    def __init__(self, feature, instance):
        self.feature = feature
        self.instance = instance
    
    def __hash__(self):
        return hash((self.feature, self.instance))
    
    def __eq__(self, other):
        return (self.feature, self.instance) == (other.feature, other.instance)
    
    def __str__(self):
        return f"<{self.feature},{self.instance}>"

class FeatureInstance:
    def __init__(self, position, instanceid):
        self.pos = position
        self.id = instanceid

class BoundingBox:
    def __init__(self, min_x, min_y, max_x, max_y):
        self.min_x = min_x
        self.min_y = min_y
        self.max_x = max_x
        self.max_y = max_y

    def intersects(self, other):
        return (self.min_x <= other.max_x and self.max_x >= other.min_x and
                self.min_y <= other.max_y and self.max_y >= other.min_y)

    @staticmethod
    def from_instances(instances, distance):
        min_x = min(inst.pos.x for inst in instances) - distance
        min_y = min(inst.pos.y for inst in instances) - distance
        max_x = max(inst.pos.x for inst in instances) + distance
        max_y = max(inst.pos.y for inst in instances) + distance
        return BoundingBox(min_x, min_y, max_x, max_y)

class ICPITree:
    def __init__(self):
        self.tree = defaultdict(lambda: defaultdict(set))

    def insert(self, instance):
        self.tree[instance.id.feature][instance.id.instance].add(instance)

    def get_neighbors(self, instance, feature, distance):
        neighbors = set()
        for other_instance in self.tree[feature].values():
            for inst in other_instance:
                if euclidean_distance(instance.pos, inst.pos) <= distance:
                    neighbors.add(inst)
        return neighbors

def euclidean_distance(pos1, pos2):
    return math.sqrt((pos1.x - pos2.x)**2 + (pos1.y - pos2.y)**2)

def generate_candidates(prev_candidates, k):
    new_candidates = set()
    for c1 in prev_candidates:
        for c2 in prev_candidates:
            if len(c1.union(c2)) == k:
                new_candidates.add(c1.union(c2))
    return new_candidates

def is_prevalent(instances, candidate, distance, min_prev):
    icpi_tree = ICPITree()
    for inst in instances:
        icpi_tree.insert(inst)

    feature_counts = defaultdict(int)
    participation_counts = defaultdict(int)

    for instance in instances:
        if instance.id.feature in candidate:
            feature_counts[instance.id.feature] += 1
            is_participating = all(
                icpi_tree.get_neighbors(instance, feature, distance)
                for feature in candidate if feature != instance.id.feature
            )
            if is_participating:
                participation_counts[instance.id.feature] += 1

    prevalences = [
        participation_counts[feature] / feature_counts[feature]
        for feature in candidate
        if feature_counts[feature] > 0
    ]

    # If prevalences is empty, it means no instances of the candidate features were found
    if not prevalences:
        return False

    return min(prevalences) >= min_prev

def BBmaxspatiotempcolloc(input_data, maxdist, minprev, minfreq, distance_metric=None, predictions=False, verbose=0, clean_trees=False):
    result = []
    all_features = set()
    for state in input_data.states.values():
        all_features.update(inst.id.feature for inst in state.instances)

    # Step 1 & 2: Find prevalent pairs and remove non-time-prevalent pairs
    prevalent_pairs = defaultdict(set)
    for time, state in input_data.states.items():
        for pair in combinations(all_features, 2):
            if is_prevalent(state.instances, set(pair), maxdist, minprev):
                prevalent_pairs[time].add(frozenset(pair))

    time_prevalent_pairs = {
        pair for pair in set.union(*prevalent_pairs.values())
        if sum(pair in pairs for pairs in prevalent_pairs.values()) >= minfreq * len(input_data.states)
    }

    # Step 3 & 4: Generate and prune candidate set
    candidates = time_prevalent_pairs
    k = 3

    while candidates:
        new_candidates = generate_candidates(candidates, k)
        prevalent_candidates = set()

        for candidate in new_candidates:
            time_prevalent_count = 0
            for time, state in input_data.states.items():
                if is_prevalent(state.instances, candidate, maxdist, minprev):
                    time_prevalent_count += 1
            
            if time_prevalent_count >= minfreq * len(input_data.states):
                prevalent_candidates.add(candidate)
                result.append(candidate)

        candidates = prevalent_candidates
        k += 1

    # Ensure maximality
    maximal_result = []
    for pattern in sorted(result, key=len, reverse=True):
        if not any(pattern.issubset(maximal) for maximal in maximal_result):
            maximal_result.append(pattern)
    
    return [frozenset(pattern) for pattern in maximal_result]

# This function signature matches what's expected in main.py
def maxspatiotempcolloc(input_data, maxdist, minprev, minfreq, predictions=False, save_trees=0, verbose=0, clean_trees=False):
    return BBmaxspatiotempcolloc(input_data, maxdist, minprev, minfreq, None, predictions, verbose, clean_trees)