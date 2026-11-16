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
def _bounded_adjacency(values, threshold):
    # bounds: L <= SMG distance <= 2L, with v = row i - row j and L = d * sum(|v - mean(v)|)
    n, dimensions = values.shape
    means = np.empty(n, dtype=np.float64)
    for i in prange(n):
        means[i] = np.sum(values[i]) / dimensions
    adjacency = np.zeros((n, n), dtype=np.bool_)
    # margin for float64 rounding
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
                # lower bound above the threshold: no edge
                if lower > threshold + margin:
                    edge = False
                # upper bound below the threshold, or rows differ only by a constant: edge
                elif 2 * lower < threshold - margin or deviation == 0:
                    edge = True
                # otherwise the exact distance
                else:
                    edge = slope_distance(values[i], values[j]) <= threshold
                adjacency[i, j] = edge
                adjacency[j, i] = edge
    return adjacency


def adjacency_matrix(values, threshold):
    # adjacency (SMG distance <= threshold) of min-max scaled rows
    return _bounded_adjacency(np.ascontiguousarray(values, dtype=np.float64), threshold)


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


def changed_degrees(first, second):
    # OTU IDs of both files, matched by ID
    names = list(dict.fromkeys([*first[0], *second[0]]))
    positions = {name: i for i, name in enumerate(names)}
    # edges in one network but not the other
    changed = np.zeros((len(names), len(names)), dtype=bool)
    for ids, adjacency in (first, second):
        indices = [positions[name] for name in ids]
        changed[np.ix_(indices, indices)] ^= adjacency
    # changed edges of each OTU
    return np.asarray(names, dtype=object), changed.sum(axis=1)


def write_changed_degrees(path, names, degrees):
    # total number of changed edges, then each OTU's count, largest first
    rows = sorted(zip(names, degrees), key=lambda row: row[1], reverse=True)
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, lineterminator='\n')
        writer.writerow(['Direct Comparison', '', ''])
        writer.writerow(['Edges', int(degrees.sum() // 2), ''])
        writer.writerows((name, int(degree), '') for name, degree in rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare two OTU networks and export one result CSV")
    parser.add_argument('-in1', required=True, help='First headerless OTU input CSV')
    parser.add_argument('-in2', required=True, help='Second headerless OTU input CSV')
    parser.add_argument('-o', '--output', required=True, help='Result CSV file')
    parser.add_argument('-pos1', type=nonnegative_float,
                        help='Fixed threshold for file 1 (default: its adaptive threshold)')
    parser.add_argument('-pos2', type=nonnegative_float,
                        help='Fixed threshold for file 2 (default: its adaptive threshold)')
    add_network_arguments(parser)
    args = parser.parse_args(argv)
    try:
        output = validate_csv_output(args.output, [args.in1, args.in2])
        networks = []
        for path, fixed in ((args.in1, args.pos1), (args.in2, args.pos2)):
            data = read_input_csv(path)
            names = data.iloc[:, 0].to_numpy()
            scaled = prepare_values(data.iloc[:, 1:].to_numpy(), args.column_normalize)
            threshold = group_threshold(scaled, args.angle, args.lam, fixed)
            networks.append((names, adjacency_matrix(scaled, threshold)))
        names, degrees = changed_degrees(*networks)
        output.parent.mkdir(parents=True, exist_ok=True)
        write_changed_degrees(output, names, degrees)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f'DAA results: {output}')


if __name__ == '__main__':
    main()
