from collections import defaultdict
import re
import pandas as pd
import numpy as np
import sys
import os
import colorsys
from openpyxl.styles import PatternFill
import argparse
import multiprocessing
from multiprocessing import shared_memory
from numba import njit
import time

@njit
def slope_distance_numba(row1, row2):
    v = row1 - row2
    v_sorted = np.sort(v)
    n = len(v_sorted)
    distance = 0
    for k in range(n):
        distance += (2 * k - n + 1) * v_sorted[k]
    return 2.0 * distance

def generate_unique_colors(num_clusters):
    colors = []
    for i in range(num_clusters):
        hue = (0.618033988749895 * i) % 1.0
        saturation = 0.5 + (0.5 * (i % 2))
        value = 0.7 + (0.3 * ((i // 2) % 2))
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        hex_color = ''.join([f'{int(x * 255):02X}' for x in rgb])
        colors.append(hex_color)
    return colors

def extract_numeric_suffix(name: str) -> int:
    if isinstance(name, bytes):
        name = name.decode('utf-8', errors='ignore')
    else:
        name = str(name)
    match = re.search(r'(\d+)$', name)
    if match:
        return int(match.group(1))
    return 0

def process_row_shared(args):
    shm_name, shape, dtype, data_global, core_id, total_cores, threshold = args
    
    existing_shm = shared_memory.SharedMemory(name=shm_name)
    dist_matrix = np.ndarray(shape, dtype=dtype, buffer=existing_shm.buf)
    
    n = shape[0]
    local_uf_edges = [] 

    for i in range(core_id, n, total_cores):
        row_i = data_global[i]
        for j in range(i + 1, n):
            dist = slope_distance_numba(row_i, data_global[j])
            
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
            
            if dist < threshold:
                local_uf_edges.append((i, j))
                
    existing_shm.close()
    return local_uf_edges

def parallel_slope_distance_shared(df_values, threshold):
    n, m = df_values.shape
    
    try:
        shm = shared_memory.SharedMemory(create=True, size=n * n * 8)
    except OSError:
        print("Error: Memory allocation failed.")
        sys.exit(1)

    shared_matrix = np.ndarray((n, n), dtype=np.float64, buffer=shm.buf)
    shared_matrix.fill(0.0)

    cpu_core_count = os.cpu_count() or 4
    
    args_list = []
    for core_id in range(cpu_core_count):
        args_list.append((
            shm.name, (n, n), np.float64, 
            df_values,
            core_id, 
            cpu_core_count, 
            threshold
        ))

    uf = UnionFind(n)
    with multiprocessing.Pool(cpu_core_count) as pool:
        all_edges_lists = pool.map(process_row_shared, args_list)

    for edges in all_edges_lists:
        for i, j in edges:
            uf.union(i, j)

    final_matrix = shared_matrix.copy()

    shm.close()
    shm.unlink()

    return final_matrix, uf

class UnionFind:
    def __init__(self, n):
        self.parent = np.arange(n, dtype=np.int32)
        self.rank = np.zeros(n, dtype=np.int8)

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        root_x, root_y = self.find(x), self.find(y)
        if root_x == root_y:
            return
        if self.rank[root_x] < self.rank[root_y]:
            root_x, root_y = root_y, root_x
        self.parent[root_y] = root_x
        if self.rank[root_x] == self.rank[root_y]:
            self.rank[root_x] += 1

def check_blank(data, file_name):
    if data.iloc[:, 0].isnull().any():
        raise ValueError(f"Blank row in '{file_name}'")

def check_duplicate(data, file_name):
    if data.iloc[:, 0].duplicated().any():
        raise ValueError(f"Duplicate object in '{file_name}'")

def min_max_scale(x):
    x_min, x_max = x.min(), x.max()
    return (x - x_min) / (x_max - x_min) if x_max != x_min else x * 0

def save_colored_network(dist_matrix, groups, threshold, output_file, names):
    colors = generate_unique_colors(len(groups))
    group_list = list(groups.items())
    group_list.sort(key=lambda x: min(extract_numeric_suffix(m) for m in x[1]))
    new_order = []
    for _, members in group_list:
        members.sort(key=lambda x: extract_numeric_suffix(x))
        new_order.extend(members)
    
    if not isinstance(dist_matrix, pd.DataFrame):
        dist_matrix = pd.DataFrame(dist_matrix, index=names, columns=names)
        
    df_dist = dist_matrix.loc[new_order, new_order]

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df_dist.to_excel(writer, sheet_name='Distance Matrix', index=True)
        worksheet = writer.book['Distance Matrix']
        
        fills = [PatternFill(start_color=c, end_color=c, fill_type='solid') for c in colors]
        
        name_to_idx = {name: i for i, name in enumerate(df_dist.index)}
        
        for group_idx, (_, members) in enumerate(groups.items()):
            if len(members) < 2: continue
            fill = fills[group_idx % len(fills)]
            
            for member in members:
                if member not in name_to_idx: continue
                r_idx = name_to_idx[member] + 2
                
                worksheet.cell(row=r_idx, column=1).fill = fill
                
                member_idx_in_matrix = name_to_idx[member]
                for other in members:
                    if other == member: continue
                    if other not in name_to_idx: continue
                    
                    other_idx_in_matrix = name_to_idx[other]
                    val = df_dist.iloc[member_idx_in_matrix, other_idx_in_matrix]
                    
                    if val < threshold:
                        c_idx = name_to_idx[other] + 2
                        worksheet.cell(row=r_idx, column=c_idx).fill = fill

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-t', '--threshold', type=float)
    args = parser.parse_args()

    try:
        col_num = len(pd.read_csv(args.input, header=None).columns)
        df = pd.read_csv(args.input, header=None, dtype={i: float for i in range(1, col_num)})
        threshold = args.threshold if args.threshold is not None else (col_num - 2)
        print(f"Using threshold: {threshold}")

        names = df.iloc[:, 0].values
        df_values = df.iloc[:, 1:]
        check_blank(df, args.input)
        check_duplicate(df, args.input)
        df_values = df_values.apply(min_max_scale, axis=1)
        
        data_matrix = df_values.to_numpy(dtype=np.float64)

        dist_matrix, uf = parallel_slope_distance_shared(data_matrix, threshold)

        n = data_matrix.shape[0]
        groups_dict = defaultdict(list)
        for i in range(n):
            root = uf.find(i)
            groups_dict[root].append(names[i])

        save_colored_network(dist_matrix, groups_dict, threshold=threshold, output_file=args.output, names=names)

        print("\nClustering Results:")
        for i, (key, group) in enumerate(groups_dict.items(), 1):
            if len(group) >= 2:
                print(f"Cluster {i}: {len(group)} members")

    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    start = time.perf_counter()
    multiprocessing.freeze_support()
    main()
    end = time.perf_counter()
    print(f"Time: {end - start:.4f}s")
