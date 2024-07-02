def manhattan_distance(point1, point2):
    return sum(abs(a - b) for a, b in zip(point1, point2))

def chebyshev_distance(point1, point2):
    return max(abs(a - b) for a, b in zip(point1, point2))

class BoundingBox:
    def __init__(self, min_coords, max_coords):
        self.min_coords = min_coords
        self.max_coords = max_coords

    def intersects(self, other_box):
        for i in range(len(self.min_coords)):
            if self.max_coords[i] < other_box.min_coords[i] or self.min_coords[i] > other_box.max_coords[i]:
                return False
        return True

    def contains_point(self, point):
        for i in range(len(point)):
            if point[i] < self.min_coords[i] or point[i] > self.max_coords[i]:
                return False
        return True
class CoLocationInstance:
    def __init__(self, points, metric='manhattan'):
        self.points = points
        self.metric = metric
        self.bounding_box = self.create_bounding_box()

    def create_bounding_box(self):
        min_coords = [min(coords) for coords in zip(*self.points)]
        max_coords = [max(coords) for coords in zip(*self.points)]
        return BoundingBox(min_coords, max_coords)

    def is_co_located(self, other_instance):
        return self.bounding_box.intersects(other_instance.bounding_box)

    def __repr__(self):
        return f"CoLocationInstance(points={self.points}, metric={self.metric})"

def generate_candidates(patterns):
    candidates = []
    for pattern in patterns:
        # Generate bounding box for the pattern
        colocation_instance = CoLocationInstance(pattern)
        candidates.append(colocation_instance)
    return candidates

class MAXMDCOPMiner:
    def __init__(self, data, metric='manhattan', distance_threshold=1.0):
        self.data = data
        self.metric = metric
        self.distance_threshold = distance_threshold

    def find_spatially_prevalent_patterns(self, ST, R):
        C2 = {}
        for Spt in ST:
            C2[tuple(map(tuple, Spt))] = self.compute_C2(Spt, R)
        return C2


    def compute_C2(self, Spt, R):
        patterns = []
        for i, point1 in enumerate(Spt):
            for j, point2 in enumerate(Spt):
                if i < j:
                    if self.metric == 'manhattan':
                        distance = manhattan_distance(point1, point2)
                    else:
                        distance = chebyshev_distance(point1, point2)
                    if distance <= R:
                        patterns.append([point1, point2])
        return patterns

    def compute_time_prevalence(self, C2, ST, mintprev):
        MDCOP_2 = {}
        for Spt, patterns in C2.items():
            MDCOP_2[Spt] = self.filter_time_prevalent(patterns, ST, mintprev)
        return MDCOP_2

    def filter_time_prevalent(self, patterns, ST, mintprev):
        return patterns

    def generate_local_candidates(self, MDCOP_2, ST):
        local_candidates = {}
        for Spt, patterns in MDCOP_2.items():
            local_candidates[Spt] = generate_candidates(patterns)
        return local_candidates

    def construct_global_candidates(self, local_candidates):
        global_candidates = set()
        for candidates in local_candidates.values():
            for candidate in candidates:
                global_candidates.add(tuple(map(tuple, candidate.points)))
        return [CoLocationInstance(list(map(list, candidate))) for candidate in global_candidates]

    def mine_maximal_MDCOPs(self, global_candidates, ST, minprev):
        maximal_MDCOPs = []
        for candidate in global_candidates:
            if self.is_time_prevalent(candidate, ST, minprev):
                maximal_MDCOPs.append(candidate)
            else:
                subsets = self.generate_subsets(candidate)
                for subset in subsets:
                    if self.is_time_prevalent(subset, ST, minprev):
                        maximal_MDCOPs.append(subset)
        return maximal_MDCOPs

    def is_time_prevalent(self, candidate, ST, minprev):
        return True

    def generate_subsets(self, candidate):
        subsets = []
        for i in range(len(candidate.points)):
            subset = candidate.points[:i] + candidate.points[i+1:]
            subsets.append(CoLocationInstance(subset))
        return subsets

    def run(self):
        ST = [self.data]
        R = self.distance_threshold
        C2 = self.find_spatially_prevalent_patterns(ST, R)
        MDCOP_2 = self.compute_time_prevalence(C2, ST, self.distance_threshold)
        local_candidates = self.generate_local_candidates(MDCOP_2, ST)
        global_candidates = self.construct_global_candidates(local_candidates)
        maximal_MDCOPs = self.mine_maximal_MDCOPs(global_candidates, ST, self.distance_threshold)
        return maximal_MDCOPs

if __name__ == '__main__':
    # Example usage
    data = [
        [1, 1, 1, 0, 0],
        [1, 2, 1, 1, 0],
        [1, 3, 1, 1, 1],
        [1, 4, 1, 0, 1],
        [1, 5, 1, 5, 5],
        [1, 6, 1, 6, 6],
        [2, 1, 1, 0, 0],
        [2, 2, 1, 1, 0],
        [2, 3, 1, 1, 1],
        [2, 4, 1, 0, 2],
        [2, 5, 1, 5, 5],
        [2, 6, 1, 6, 6],
        [3, 1, 1, 0, 0],
        [3, 2, 1, 1, 0],
        [3, 3, 1, 1, 1],
        [3, 4, 1, 2, 0]
    ]

    miner = MAXMDCOPMiner(data, metric='manhattan', distance_threshold=1.0)
    co_location_patterns = miner.run()
    print(co_location_patterns)