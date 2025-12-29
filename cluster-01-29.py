import argparse
import sys
import time
import re
import random
from collections import defaultdict

import pandas as pd
import numpy as np
from openpyxl import Workbook
from openpyxl.styles import PatternFill

from joblib import Parallel, delayed
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from numba import njit


@njit
def slope_distance(row1, row2):
    v = row1 - row2
    v_sorted = np.sort(v)
    n = len(v_sorted)
    distance = 0.0
    for k in range(n):
        distance += (2 * k - n + 1) * v_sorted[k]
    return 2.0 * distance


def extract_numeric_suffix(name):
    try:
        return int(re.findall(r'\d+', str(name))[-1])
    except (IndexError, ValueError, TypeError):
        return 0


def calculate_row_edges(i, data_matrix, threshold):
    n = data_matrix.shape[0]
    row_i = data_matrix[i]
    edges = []

    for j in range(i + 1, n):
        dist = slope_distance(row_i, data_matrix[j])
        if dist < threshold:
            edges.append((i, j, dist))

    return edges


def get_clusters(data_matrix, threshold):
    n = data_matrix.shape[0]
    print(f"Calculating distances for {n} items using parallel processing...")

    results = Parallel(n_jobs=-1, backend='loky')(
        delayed(calculate_row_edges)(i, data_matrix, threshold)
        for i in range(n)
    )

    all_edges = [edge for row_res in results for edge in row_res]

    rows = [e[0] for e in all_edges]
    cols = [e[1] for e in all_edges]
    data = [1] * len(all_edges)

    graph = csr_matrix((data, (rows, cols)), shape=(n, n))

    n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)

    dist_matrix_out = np.zeros((n, n))
    for i, j, d in all_edges:
        dist_matrix_out[i, j] = d
        dist_matrix_out[j, i] = d

    return dist_matrix_out, labels


def check_data_integrity(df, filename):
    if df.iloc[:, 0].isnull().any():
        raise ValueError(f"Found blank rows in '{filename}'")
    if df.iloc[:, 0].duplicated().any():
        raise ValueError(f"Found duplicate identifiers in '{filename}'")


def min_max_scale(x):
    diff = x.max() - x.min()
    if diff == 0:
        return x * 0
    return (x - x.min()) / diff


def get_color_list(n):
    presets = [
        "FFB6C1", "87CEFA", "90EE90", "F0E68C", "D3D3D3", "FFA07A",
        "20B2AA", "9370DB", "FF6347", "FFD700", "B0C4DE", "FFDEAD"
    ]

    colors = presets[:]
    while len(colors) < n:
        colors.append("".join([random.choice('89ABCDEF') for _ in range(6)]))
    return colors


def save_excel(dist_matrix, groups, threshold, output_file, names):
    print("Saving results to Excel...")

    group_items = sorted(groups.items(), key=lambda x: min(extract_numeric_suffix(n) for n in x[1]))

    new_order_indices = []
    name_to_idx = {name: i for i, name in enumerate(names)}

    for _, members in group_items:
        members.sort(key=extract_numeric_suffix)
        for m in members:
            if m in name_to_idx:
                new_order_indices.append(name_to_idx[m])

    df_dist = pd.DataFrame(dist_matrix, index=names, columns=names)
    df_dist = df_dist.iloc[new_order_indices, new_order_indices]

    colors = get_color_list(len(groups))

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df_dist.to_excel(writer, sheet_name='Clustering', index=True)
        ws = writer.book['Clustering']

        fills = [PatternFill(start_color=c, end_color=c, fill_type='solid') for c in colors]

        excel_name_to_idx = {name: i for i, name in enumerate(df_dist.index)}

        for group_idx, (_, members) in enumerate(group_items):
            if len(members) < 2:
                continue

            fill = fills[group_idx % len(fills)]

            for member in members:
                if member not in excel_name_to_idx: continue

                r = excel_name_to_idx[member] + 2

                ws.cell(row=r, column=1).fill = fill

                member_idx_orig = name_to_idx[member]

                for other in members:
                    if member == other: continue
                    if other not in excel_name_to_idx: continue

                    other_idx_orig = name_to_idx[other]
                    val = dist_matrix[member_idx_orig, other_idx_orig]

                    if val < threshold:
                        c = excel_name_to_idx[other] + 2
                        ws.cell(row=r, column=c).fill = fill


def main():
    parser = argparse.ArgumentParser(description="Cluster Analysis Tool")
    parser.add_argument('-i', '--input', required=True, help="Input CSV file")
    parser.add_argument('-o', '--output', required=True, help="Output Excel file")
    parser.add_argument('-t', '--threshold', type=float, help="Clustering threshold")
    args = parser.parse_args()

    try:
        preview = pd.read_csv(args.input, header=None, nrows=1)
        col_num = len(preview.columns)

        df = pd.read_csv(args.input, header=None)

        check_data_integrity(df, args.input)

        names = df.iloc[:, 0].values
        df_values = df.iloc[:, 1:].astype(float).apply(min_max_scale, axis=1)
        data_matrix = df_values.to_numpy()

        threshold = args.threshold if args.threshold is not None else (col_num - 2)
        print(f"Using threshold: {threshold}")

        dist_matrix, labels = get_clusters(data_matrix, threshold)

        groups = defaultdict(list)
        for idx, label in enumerate(labels):
            groups[label].append(names[idx])

        save_excel(dist_matrix, groups, threshold, args.output, names)

        print("\nSummary:")
        valid_groups = [g for g in groups.values() if len(g) >= 2]
        print(f"Found {len(valid_groups)} clusters with >= 2 members.")
        for i, g in enumerate(valid_groups, 1):
            print(f"  Cluster {i}: {len(g)} members")

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    start_time = time.perf_counter()
    main()
    print(f"Total time: {time.perf_counter() - start_time:.2f}s")