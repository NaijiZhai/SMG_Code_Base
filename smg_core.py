"""SMG distances, adaptive thresholds and undirected OTU networks."""
from dataclasses import dataclass
import math

import numpy as np
from numba import njit, prange


@njit(cache=True, fastmath=True, nogil=True)
def slope_distance(row1, row2):
    """Sum pairwise slope differences using sorted differences in O(n log n)."""
    differences = np.sort(row1 - row2)
    n = len(differences)
    total = 0.0
    for k in range(n):
        total += (2 * k - n + 1) * differences[k]
    return 2.0 * total


@njit(cache=True, parallel=True)
def _distance_matrix(values):
    n = len(values)
    result = np.zeros((n, n), dtype=np.float64)
    for task in prange((n + 1) // 2):
        index = np.int64(task)
        # Complementary rows contain n - 1 pairs in total.
        for side in range(2):
            i = index if side == 0 else n - 1 - index
            if side == 1 and i == index:
                continue
            for j in range(i + 1, n):
                value = slope_distance(values[i], values[j])
                result[i, j] = value
                result[j, i] = value
    return result


def distance_matrix(values):
    """Symmetric row-to-row distances using the reference sorting formula."""
    return _distance_matrix(np.ascontiguousarray(values, dtype=np.float64))


@njit(cache=True, parallel=True)
def _bounded_adjacency(values, threshold):
    """L <= SMG distance <= 2L, where L = d * sum(abs(v - mean(v)))."""
    n, dimensions = values.shape
    means = np.empty(n, dtype=np.float64)
    for i in prange(n):
        means[i] = np.sum(values[i]) / dimensions
    adjacency = np.zeros((n, n), dtype=np.bool_)
    # This conservative guard covers accumulated float64 roundoff for
    # min-max scaled vectors. Ambiguous pairs use the reference kernel.
    margin = 128 * np.finfo(np.float64).eps * dimensions**3
    for task in prange((n + 1) // 2):
        index = np.int64(task)
        for side in range(2):
            i = index if side == 0 else n - 1 - index
            if side == 1 and i == index:
                continue
            for j in range(i + 1, n):
                center = means[i] - means[j]
                deviation = 0.0
                for k in range(dimensions):
                    deviation += abs(values[i, k] - values[j, k] - center)
                lower = dimensions * deviation
                if lower > threshold + margin:
                    edge = False
                elif 2 * lower < threshold - margin or deviation == 0:
                    edge = True
                else:
                    edge = slope_distance(values[i], values[j]) <= threshold
                adjacency[i, j] = edge
                adjacency[j, i] = edge
    return adjacency


def adjacency_matrix(values, threshold):
    """Exact threshold decisions for min-max scaled rows, without a distance matrix."""
    return _bounded_adjacency(np.ascontiguousarray(values, dtype=np.float64), threshold)


def prepare_values(values, column_normalize=False):
    """Optionally normalize sample columns, then min-max scale each OTU row."""
    values = np.asarray(values, dtype=np.float64)
    if column_normalize:
        totals = values.sum(axis=0, keepdims=True)
        values = values / np.where(totals == 0, 1.0, totals)
    low = values.min(axis=1, keepdims=True)
    span = values.max(axis=1, keepdims=True) - low
    return np.ascontiguousarray((values - low) / np.where(span == 0, 1.0, span))


def arc_values(x, angle):
    """Reference arc at x in [0, 1], for a validated angle in (0, 90] degrees."""
    t = np.tan(np.deg2rad(angle) / 2)
    return 2*x*(1-t+t*x) / (1+t+np.sqrt((1-t)**2+4*t*(1-x)*(1+t*x)))


@dataclass(frozen=True)
class Threshold:
    value: float
    adaptive: float | None
    line_distance: float
    standard: np.ndarray
    deviated: np.ndarray
    angle: float
    multiplier: float
    kind: str


def group_threshold(values, angle=90.0, multiplier=1.0, fixed=None):
    """Resolve a group's threshold and retain its reference curves for reporting."""
    total = np.sort(values, axis=1).sum(axis=0)
    span = total.max() - total.min()
    if span == 0:
        if fixed is None:
            raise ValueError('Every OTU is constant; supply a fixed threshold because there is no standard OTU')
        standard = np.full(total.size, np.nan)
        deviated = standard.copy()
        adaptive = None
    else:
        standard = (total - total.min()) / span
        deviated = arc_values(standard, angle)
        adaptive = slope_distance(standard, deviated)
    value = float(adaptive if fixed is None else fixed) * multiplier
    if not math.isfinite(value) or value < 0:
        raise ValueError('the multiplied threshold must be finite and non-negative')
    line = np.linspace(0.0, 1.0, total.size)
    return Threshold(value, adaptive, slope_distance(line, arc_values(line, angle)),
                     standard, deviated, angle, multiplier, 'adaptive' if fixed is None else 'custom')


@dataclass(frozen=True)
class Network:
    names: np.ndarray
    values: np.ndarray
    distances: np.ndarray | None
    adjacency: np.ndarray
    threshold: Threshold


def build_network(names, values, *, angle=90.0, multiplier=1.0, fixed=None, column_normalize=False,
                  include_distances=False):
    scaled = prepare_values(values, column_normalize)
    threshold = group_threshold(scaled, angle, multiplier, fixed)
    if include_distances:
        distances = distance_matrix(scaled)
        adjacency = distances <= threshold.value
        np.fill_diagonal(adjacency, False)
    else:
        distances = None
        adjacency = adjacency_matrix(scaled, threshold.value)
    return Network(np.asarray(names, dtype=object), scaled, distances, adjacency, threshold)


def changed_degrees(first, second):
    """Align OTU IDs and count incident edges in the networks' symmetric difference."""
    names = list(dict.fromkeys([*first.names, *second.names]))
    positions = {name: i for i, name in enumerate(names)}
    changed = np.zeros((len(names), len(names)), dtype=bool)
    for network in (first, second):
        indices = [positions[name] for name in network.names]
        changed[np.ix_(indices, indices)] ^= network.adjacency
    return np.asarray(names, dtype=object), changed.sum(axis=1)


def connected_components(adjacency):
    """Components and their members in input order, including singleton OTUs."""
    seen = np.zeros(len(adjacency), dtype=bool)
    groups = []
    for start in range(len(adjacency)):
        if seen[start]:
            continue
        seen[start] = True
        pending, members = [start], []
        while pending:
            node = pending.pop()
            members.append(node)
            neighbors = np.flatnonzero(adjacency[node] & ~seen)
            seen[neighbors] = True
            pending.extend(neighbors)
        groups.append(sorted(members))
    return groups
