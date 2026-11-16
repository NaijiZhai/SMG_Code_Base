import argparse
import csv
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


def write_distance_matrix(path, names, distances, order=None):
    # labeled rows and columns in the same order, full float precision
    if order is None:
        order = np.arange(len(names))
    with open(path, 'w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['sample_id', *names[order]])
        for i in order:
            writer.writerow([names[i], *distances[i, order]])


def main(argv=None):
    parser = argparse.ArgumentParser(description='Compute sample-to-sample SMG distances and export one CSV')
    parser.add_argument('-i', '--input', required=True, help='Headerless CSV with one sample per row')
    parser.add_argument('-o', '--output', required=True, help='Sample distance matrix CSV')
    args = parser.parse_args(argv)
    try:
        output = validate_csv_output(args.output, [args.input])
        data = read_input_csv(args.input)
        names = data.iloc[:, 0].to_numpy()
        distances = distance_matrix(data.iloc[:, 1:].to_numpy())
        output.parent.mkdir(parents=True, exist_ok=True)
        write_distance_matrix(output, names, distances)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f'Sample distance matrix: {output}')


if __name__ == '__main__':
    main()
