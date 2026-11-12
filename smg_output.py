"""Result file writers for the SMG command-line tools."""
import colorsys
import csv
import json
from itertools import chain
from pathlib import Path
import re

import numpy as np
import pandas as pd
from openpyxl.styles import PatternFill


def adaptive_text(threshold):
    return f'{threshold.adaptive:.6f}' if threshold.adaptive is not None else 'unavailable (all OTUs constant)'


def threshold_context(threshold):
    return (f'lambda {threshold.multiplier:g}; adaptive {adaptive_text(threshold)}; '
            f'straight line with the {threshold.angle:g} degree arc {threshold.line_distance:.6f}')


def threshold_curve(threshold):
    n = len(threshold.standard)
    return pd.DataFrame({'rank': np.arange(1, n + 1), 'line': np.linspace(0.0, 1.0, n),
                         'standard_otu': threshold.standard, 'deviated_otu': threshold.deviated})


def write_thresholds(output, thresholds):
    output = Path(output)
    path = output.with_name(output.stem + '_standard_otu.csv')
    curves = [threshold_curve(threshold).assign(file=i, threshold=threshold.value)
              for i, threshold in enumerate(thresholds, 1)]
    pd.concat(curves)[['file', 'rank', 'line', 'standard_otu', 'deviated_otu', 'threshold']].to_csv(path, index=False)
    return path


def write_changed_degrees(path, names, degrees):
    rows = sorted(zip(names, degrees), key=lambda row: row[1], reverse=True)
    records = chain([('Edges', int(degrees.sum() // 2), '')],
                    ((name, int(degree), '') for name, degree in rows))
    _write_csv(path, ['Direct Comparison', '', ''], records, lineterminator='\n')


def _write_csv(path, header, rows, *, lineterminator='\r\n'):
    """Write a header and stream rows, preserving each output's newline convention."""
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, lineterminator=lineterminator)
        writer.writerow(header)
        writer.writerows(rows)


def write_distances(path, names, distances):
    pd.DataFrame(distances, index=names, columns=names).to_excel(
        path, sheet_name='Distance Matrix', engine='openpyxl', index=True)


def _write_ids(directory, names):
    _write_csv(directory / 'ids.csv', ['id'], ((name,) for name in names))


def write_binary_distances(path, names, distances):
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    np.save(directory / 'distances.npy', distances, allow_pickle=False)
    _write_ids(directory, names)


def _threshold_info(network):
    threshold = network.threshold
    return {
        'Threshold used': threshold.value, 'Kind': threshold.kind, 'Lambda': threshold.multiplier,
        'Adaptive threshold': threshold.adaptive, 'Angle': threshold.angle,
        'Straight line threshold': threshold.line_distance, 'Samples': network.values.shape[1],
        'OTUs': len(network.names), 'Non-constant OTUs': int((np.ptp(network.values, axis=1) > 0).sum()),
    }


def write_binary_clusters(path, network, groups):
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    names = network.names
    edges = ((name, names[j + i + 1]) for i, name in enumerate(names)
             for j in np.flatnonzero(network.adjacency[i, i + 1:]))
    _write_csv(directory / 'edges.csv', ['source', 'target'], edges)
    labels = np.empty(len(names), dtype=int)
    for label, members in enumerate(groups, 1):
        labels[members] = label
    _write_csv(directory / 'clusters.csv', ['otu', 'cluster'], zip(names, labels))
    (directory / 'threshold.json').write_text(json.dumps(_threshold_info(network), indent=2) + '\n', encoding='utf-8')
    threshold_curve(network.threshold).to_csv(directory / 'standard_otu.csv', index=False)
    if network.distances is not None:
        np.save(directory / 'distances.npy', network.distances, allow_pickle=False)
        _write_ids(directory, names)
    else:
        (directory / 'distances.npy').unlink(missing_ok=True)
        (directory / 'ids.csv').unlink(missing_ok=True)


def _numeric_suffix(name):
    match = re.search(r'(\d+)$', str(name))
    return int(match.group(1)) if match else 0


def _cluster_fill(index):
    rgb = colorsys.hsv_to_rgb((0.618033988749895 * index) % 1.0,
                             0.5 + 0.5 * (index % 2), 0.7 + 0.3 * ((index // 2) % 2))
    color = ''.join(f'{int(value * 255):02X}' for value in rgb)
    return PatternFill(start_color=color, end_color=color, fill_type='solid')


def write_clusters(path, network, groups):
    if network.distances is None:
        raise ValueError('Excel clustering output requires the full distance matrix')
    names, threshold = network.names, network.threshold
    suffix = lambda i: _numeric_suffix(names[i])
    ordered_groups = sorted(groups, key=lambda members: min(suffix(i) for i in members))
    order = [i for members in ordered_groups for i in sorted(members, key=suffix)]
    distances = pd.DataFrame(network.distances[np.ix_(order, order)], index=names[order], columns=names[order])
    info = _threshold_info(network)
    if info['Adaptive threshold'] is None:
        info['Adaptive threshold'] = adaptive_text(threshold)
    info = pd.Series(info, name='Value')
    curve = threshold_curve(threshold).rename(columns={
        'rank': 'Rank', 'line': 'Line', 'standard_otu': 'Standard OTU', 'deviated_otu': 'Deviated OTU',
    })
    positions = np.empty(len(names), dtype=int)
    positions[order] = np.arange(2, len(names) + 2)
    with pd.ExcelWriter(path, engine='openpyxl') as writer:
        distances.to_excel(writer, sheet_name='Distance Matrix', index=True)
        info.to_excel(writer, sheet_name='Threshold', index=True)
        curve.to_excel(writer, sheet_name='Standard OTU', index=False)
        sheet = writer.book['Distance Matrix']
        for color, members in enumerate(group for group in groups if len(group) >= 2):
            fill = _cluster_fill(color)
            for member in members:
                row = int(positions[member])
                sheet.cell(row=row, column=1).fill = fill
                for other in np.flatnonzero(network.adjacency[member]):
                    sheet.cell(row=row, column=int(positions[other])).fill = fill
