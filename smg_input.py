"""Shared validation for headerless SMG inputs and network parameters."""
import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


def read_input_csv(path):
    """Read string IDs and finite, non-negative abundances without inferring IDs."""
    try:
        cells = pd.read_csv(path, header=None, dtype=str, encoding='utf-8-sig',
                            keep_default_na=False, na_filter=False, skip_blank_lines=False)
    except pd.errors.EmptyDataError as exc:
        raise ValueError(f'{path}: input is empty') from exc
    if cells.empty or cells.shape[1] < 2:
        raise ValueError(f'{path}: each row needs an ID and at least one abundance value')

    ids = cells.iloc[:, 0].str.strip()
    blank = np.flatnonzero(ids.eq(''))
    if blank.size:
        raise ValueError(f'{path}: blank ID at row {blank[0] + 1}, column 1')
    duplicates = np.flatnonzero(ids.duplicated())
    if duplicates.size:
        row = duplicates[0]
        raise ValueError(f'{path}: duplicate ID {ids.iloc[row]!r} at row {row + 1}')

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
    """Options shared by differential network analysis and clustering."""
    parser.add_argument('-a', '--angle', type=angle_degrees, default=90.0,
                        help='Arc angle in degrees (0 < angle <= 90)')
    parser.add_argument('-l', '--lam', type=nonnegative_float, default=1.0,
                        help='Finite non-negative multiplier for network thresholds')
    parser.add_argument('--column-normalize', action='store_true',
                        help='Normalize each sample column to relative abundance before row min-max')


def add_output_arguments(parser, clustering=False, default_format='binary'):
    parser.add_argument('-o', '--output', required=True, help='Output directory (binary) or .xlsx file (excel)')
    parser.add_argument('--format', choices=('binary', 'excel'), default=default_format,
                        help=f'Output format: binary result directory or Excel workbook (default: {default_format})')
    if clustering:
        parser.add_argument('--distances', action='store_true',
                            help='Also save a full distances.npy matrix in binary mode; Excel always includes distances')


def validate_output(path, output_format):
    path = Path(path)
    if output_format == 'binary':
        if path.suffix.lower() in ('.xlsx', '.xls', '.npy', '.npz'):
            raise ValueError('binary output needs a directory path; use --format excel for an .xlsx file')
        if path.exists() and not path.is_dir():
            raise ValueError(f'{path}: binary output must be a directory')
    elif path.suffix.lower() != '.xlsx' or path.is_dir():
        raise ValueError('Excel output must be an .xlsx file; use --format binary for a result directory')
