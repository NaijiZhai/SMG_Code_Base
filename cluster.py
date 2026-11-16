import argparse
import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit, prange


def read_input_csv(path):
    # read every cell as text, no header row
    try:
        cells = pd.read_csv(path, header=None, dtype=str, encoding='utf-8-sig',
                            keep_default_na=False, na_filter=False, skip_blank_lines=False)
    except pd.errors.EmptyDataError as exc:
        raise ValueError(f'{path}: input is empty') from exc
    if cells.empty or cells.shape[1] < 2:
        raise ValueError(f'{path}: each row needs an ID and at least one abundance value')

    # first column: IDs, none blank, no duplicates
    ids = cells.iloc[:, 0].str.strip()
    blank = np.flatnonzero(ids.eq(''))
    if blank.size:
        raise ValueError(f'{path}: blank ID at row {blank[0] + 1}, column 1')
    duplicates = np.flatnonzero(ids.duplicated())
    if duplicates.size:
        row = duplicates[0]
        raise ValueError(f'{path}: duplicate ID {ids.iloc[row]!r} at row {row + 1}')

    # other columns: abundances, numbers >= 0
    values = cells.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').astype(np.float64)
    bad = np.argwhere(~np.isfinite(values.to_numpy()) | (values.to_numpy() < 0))
    if len(bad):
        row, col = bad[0]
        raise ValueError(f'{path}: invalid abundance at row {row + 1}, column {col + 2} '
                         f'(ID {ids.iloc[row]!r}); values must be finite and non-negative')
    return pd.concat([ids, values], axis=1)


def angle_degrees(value):
    angle = float(value)
    if not math.isfinite(angle) or not 0 < angle <= 90:
        raise argparse.ArgumentTypeError('angle must be finite and satisfy 0 < angle <= 90')
    return angle


def nonnegative_float(value):
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError('value must be finite and non-negative')
    return number


def add_network_arguments(parser):
    # arc angle, threshold multiplier and column normalization options
    parser.add_argument('-a', '--angle', type=angle_degrees, default=90.0,
                        help='Arc angle in degrees (0 < angle <= 90)')
    parser.add_argument('-l', '--lam', type=nonnegative_float, default=1.0,
                        help='Finite non-negative multiplier for network thresholds')
    parser.add_argument('--column-normalize', action='store_true',
                        help='Normalize each sample column to relative abundance before row min-max')


def validate_csv_output(path, inputs):
    # output must be a .csv file and must not be one of the input files
    path = Path(path)
    if path.suffix.lower() != '.csv' or path.is_dir():
        raise ValueError(f'{path}: output must be a .csv file')
    if path.resolve() in {Path(source).resolve() for source in inputs}:
        raise ValueError(f'{path}: output must be different from the input files')
    return path


@njit(cache=True, fastmath=True, nogil=True)
def slope_distance(row1, row2):
    # differences between the two rows, sorted
    differences = np.sort(row1 - row2)
    n = len(differences)
    # sum of |d_k - d_l| over all pairs k < l, from the sorted order in O(n log n)
    total = 0.0
    for k in range(n):
        total += (2 * k - n + 1) * differences[k]
    # count both orders of each pair
    return 2.0 * total


@njit(cache=True, parallel=True)
def _distance_matrix(values):
    n = len(values)
    result = np.zeros((n, n), dtype=np.float64)
    for task in prange((n + 1) // 2):
        index = np.int64(task)
        # pair row i with row n - 1 - i so every task has n - 1 pairs
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
    # symmetric matrix of SMG distances between rows
    return _distance_matrix(np.ascontiguousarray(values, dtype=np.float64))


def prepare_values(values, column_normalize=False):
    values = np.asarray(values, dtype=np.float64)
    # optional: each sample column to relative abundance
    if column_normalize:
        totals = values.sum(axis=0, keepdims=True)
        values = values / np.where(totals == 0, 1.0, totals)
    # min-max scale each OTU row
    low = values.min(axis=1, keepdims=True)
    span = values.max(axis=1, keepdims=True) - low
    return np.ascontiguousarray((values - low) / np.where(span == 0, 1.0, span))


def arc_values(x, angle):
    # deviated OTU: arc from (0, 0) to (1, 1) with central angle `angle`, value at x in [0, 1]
    t = np.tan(np.deg2rad(angle) / 2)
    return 2*x*(1-t+t*x) / (1+t+np.sqrt((1-t)**2+4*t*(1-x)*(1+t*x)))


def group_threshold(values, angle=90.0, multiplier=1.0, fixed=None):
    # adaptive threshold unless a fixed one is given
    if fixed is None:
        # standard OTU: sort each row, sum by rank, min-max scale
        total = np.sort(values, axis=1).sum(axis=0)
        span = total.max() - total.min()
        if span == 0:
            raise ValueError('Every OTU is constant; supply a fixed threshold because there is no standard OTU')
        standard = (total - total.min()) / span
        # threshold = distance between the standard OTU and its deviated OTU
        fixed = slope_distance(standard, arc_values(standard, angle))
    # times the multiplier
    value = float(fixed) * multiplier
    if not math.isfinite(value) or value < 0:
        raise ValueError('the multiplied threshold must be finite and non-negative')
    return value


def connected_components(adjacency):
    # connected components of the network, single OTUs included
    seen = np.zeros(len(adjacency), dtype=bool)
    groups = []
    for start in range(len(adjacency)):
        if seen[start]:
            continue
        seen[start] = True
        # depth-first search from each OTU not yet seen
        pending, members = [start], []
        while pending:
            node = pending.pop()
            members.append(node)
            neighbors = np.flatnonzero(adjacency[node] & ~seen)
            seen[neighbors] = True
            pending.extend(neighbors)
        groups.append(sorted(members))
    return groups


def write_distance_matrix(path, names, distances, order=None):
    # labeled rows and columns in the same order, full float precision
    if order is None:
        order = np.arange(len(names))
    with open(path, 'w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['otu_id', *names[order]])
        for i in order:
            writer.writerow([names[i], *distances[i, order]])


def order_components(groups, names):
    # largest components first, ties and members ordered by OTU ID
    ordered = [sorted(members, key=lambda i: str(names[i])) for members in groups]
    return sorted(ordered, key=lambda members: (-len(members), str(names[members[0]])))


def write_cluster_members(path, names, groups):
    # one row per OTU, only components with at least two OTUs
    with open(path, 'w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['Input OTUs not listed below are isolated nodes in the network.', '', ''])
        writer.writerow(['cluster_id', 'cluster_size', 'otu_id'])
        for cluster_id, members in enumerate(groups, 1):
            if len(members) >= 2:
                writer.writerows((cluster_id, len(members), names[i]) for i in members)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build an OTU network and export two CSV files")
    parser.add_argument('-i', '--input', required=True, help='Headerless CSV with one OTU per row')
    parser.add_argument('-o', '--output', required=True,
                        help='Distance matrix CSV; the member list uses the same stem plus _clusters.csv')
    parser.add_argument('-t', '--threshold', type=nonnegative_float,
                        help="Fixed threshold (default: the group's adaptive threshold)")
    add_network_arguments(parser)
    args = parser.parse_args(argv)
    try:
        matrix_path = validate_csv_output(args.output, [args.input])
        members_path = validate_csv_output(
            matrix_path.with_name(matrix_path.stem + '_clusters.csv'), [args.input])
        data = read_input_csv(args.input)
        names = data.iloc[:, 0].to_numpy()
        scaled = prepare_values(data.iloc[:, 1:].to_numpy(), args.column_normalize)
        threshold = group_threshold(scaled, args.angle, args.lam, args.threshold)
        distances = distance_matrix(scaled)
        adjacency = distances <= threshold
        np.fill_diagonal(adjacency, False)
        groups = order_components(connected_components(adjacency), names)
        order = np.asarray([i for members in groups for i in members], dtype=np.int64)
        matrix_path.parent.mkdir(parents=True, exist_ok=True)
        write_distance_matrix(matrix_path, names, distances, order)
        write_cluster_members(members_path, names, groups)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f'Distance matrix: {matrix_path}')
    print(f'Cluster members: {members_path}')


if __name__ == '__main__':
    main()
