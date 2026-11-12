
import argparse
import csv
import os
import re
import sys
import unicodedata

import numpy as np
import pandas as pd


def read_cells(path, sheet=None):
    """Whole table as strings/objects, header row included (row 0), nothing parsed or renamed."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xlsm", ".xls"):
        return pd.read_excel(path, sheet_name=sheet if sheet is not None else 0, header=None, dtype=object)
    sep = "\t" if ext in (".tsv", ".txt", ".tab") else ","
    return pd.read_csv(path, sep=sep, header=None, dtype=str, encoding="utf-8-sig",
                       keep_default_na=False, na_filter=False)


def as_id(v):
    """Cell value -> ID string ('' for an empty cell); Excel numbers such as 114783.0 -> '114783'."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v).strip()


def load_groups(path, sheet=None):
    """[(sample, group), ...] in file order."""
    cells = read_cells(path, sheet)
    header = [as_id(v).lower() for v in cells.iloc[0]]
    if "sample" not in header or "group" not in header:
        sys.exit(f"{path}: the first row must name the columns 'sample' and 'group' (found {header})")
    cs, cg = header.index("sample"), header.index("group")
    pairs = []
    for r in range(1, len(cells)):
        s, g = as_id(cells.iat[r, cs]), as_id(cells.iat[r, cg])
        if not s and not g:
            continue
        if not s or not g:
            sys.exit(f"{path}, row {r + 1}: sample or group is empty")
        pairs.append((s, g))
    seen = {}
    for s, g in pairs:
        if s in seen:
            sys.exit(f"{path}: sample {s} is listed more than once ({seen[s]}, {g})")
        seen[s] = g
    if not pairs:
        sys.exit(f"{path}: no samples listed")
    return pairs


def load_counts(path, samples, sheet=None):
    """Feature x sample DataFrame (float) with the requested samples, in the requested order."""
    cells = read_cells(path, sheet)
    header = [as_id(v) for v in cells.iloc[0, 1:]]
    body = cells.iloc[1:]
    ids = [as_id(v) for v in body.iloc[:, 0]]
    cols = []
    for s in samples:
        hits = [k for k, h in enumerate(header) if h == s]
        if len(hits) != 1:
            sys.exit(f"{path}: sample {s} is {'missing from' if not hits else 'repeated in'} the header row")
        cols.append(hits[0] + 1)
    sub = body.iloc[:, cols]
    vals = sub.apply(lambda c: pd.to_numeric(c, errors="coerce")).to_numpy(dtype=float)
    # drop blank lines (no ID, no values), e.g. trailing rows left by Excel
    blank = np.array([i == "" and all(as_id(v) == "" for v in sub.iloc[k]) for k, i in enumerate(ids)], dtype=bool)
    ids = [i for i, b in zip(ids, blank) if not b]
    vals = vals[~blank]
    if any(i == "" for i in ids):
        sys.exit(f"{path}: a row with values has no feature ID")
    dup = pd.Index(ids)[pd.Index(ids).duplicated()].unique()
    if len(dup):
        sys.exit(f"{path}: repeated feature IDs, e.g. {list(dup[:5])}")
    bad = np.argwhere(~np.isfinite(vals) | (vals < 0))
    if len(bad):
        r, c = bad[0]
        sys.exit(f"{path}: {len(bad)} empty, non-numeric or negative cells among the used samples, "
                 f"first at feature {ids[r]}, sample {samples[c]}")
    return pd.DataFrame(vals, index=ids, columns=samples)


def clean(counts, groups, min_mean=1e-5, min_prev=0.2):
    """counts: feature x sample (all features of the raw table); groups: {group: [samples]}.
    Returns (RA of the kept features in output order, per-feature report)."""
    totals = counts.sum(axis=0)
    empty = list(totals.index[totals <= 0])
    if empty:
        sys.exit(f"samples with a total count of 0: {empty}")
    ra = counts / totals
    means = pd.DataFrame({g: ra[s].mean(axis=1) for g, s in groups.items()})
    prevs = pd.DataFrame({g: (ra[s] > 0).sum(axis=1) / len(s) for g, s in groups.items()})
    mean_max, prev_max = means.max(axis=1), prevs.max(axis=1)
    keep = (mean_max >= min_mean) & (prev_max >= min_prev)
    order = mean_max[keep].sort_values(ascending=False, kind="mergesort").index
    report = pd.concat([means.add_prefix("mean_"), prevs.add_prefix("prevalence_")], axis=1)
    report["mean_max"], report["prevalence_max"] = mean_max, prev_max
    report["pass_mean"], report["pass_prevalence"], report["kept"] = mean_max >= min_mean, prev_max >= min_prev, keep
    report["output_rank"] = pd.Series(np.arange(1, len(order) + 1), index=order).astype("Int64")
    return ra.loc[order], report


def write_input(path, ids, matrix):
    """Header-less CSV in the format of the authors' inputs (UTF-8 with BOM, CRLF); full-precision values."""
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, lineterminator="\r\n")
        for name, row in zip(ids, matrix):
            w.writerow([name] + ["0" if v == 0 else repr(float(v)) for v in row])


def safe_name(s):
    return re.sub(r'[\\/:*?"<>|\s]+', "_", s)


def main():
    ap = argparse.ArgumentParser(description="SMG data cleaning (relative abundance + mean/prevalence filter); "
                                             "writes one input file per group.")
    ap.add_argument("-c", "--counts", required=True, help="raw count table (.csv, .tsv/.txt or .xlsx)")
    ap.add_argument("--sheet", help="sheet name of an .xlsx count table (default: first sheet)")
    ap.add_argument("-g", "--groups", required=True, help="table with columns sample, group")
    ap.add_argument("--groups-sheet", help="sheet name of an .xlsx group table (default: first sheet)")
    ap.add_argument("-o", "--out-dir", required=True, help="output directory (created if missing)")
    ap.add_argument("--layout", choices=["features", "samples"], default="features",
                    help="rows of the output: features (DAA, cluster; default) or samples (BDA)")
    ap.add_argument("--prefix", default="input_", help="output file name = prefix + group + .csv (default input_)")
    ap.add_argument("--combined", metavar="NAME", help="write all groups into one file NAME.csv instead")
    ap.add_argument("--min-mean", type=float, default=1e-5,
                    help="minimum of the largest group mean relative abundance (default 0.00001 = 0.001%%)")
    ap.add_argument("--min-prev", type=float, default=0.2,
                    help="minimum of the largest group prevalence, share of samples > 0 (default 0.2)")
    args = ap.parse_args()

    pairs = load_groups(args.groups, args.groups_sheet)
    groups = {}
    for s, g in pairs:
        groups.setdefault(g, []).append(s)
    samples = [s for s, _ in pairs]
    counts = load_counts(args.counts, samples, args.sheet)
    kept, report = clean(counts, groups, args.min_mean, args.min_prev)

    n_all, n_kept = len(report), len(kept)
    print(f"raw table: {n_all} features; {len(samples)} samples in {len(groups)} groups: "
          + ", ".join(f"{g} {len(s)}" for g, s in groups.items()))
    print(f"rule: largest group mean RA >= {args.min_mean:g} and largest group prevalence >= {args.min_prev:g}")
    print(f"kept {n_kept} of {n_all} features ({100 * n_kept / n_all:.1f}%); removed: "
          f"{int((~report.pass_mean).sum())} below the mean, {int((~report.pass_prevalence).sum())} below the "
          f"prevalence (a feature can fail both); all-zero in these samples: {int((report.mean_max == 0).sum())}")
    if n_kept == 0:
        sys.exit("no feature passed the filter; nothing written")
    if args.combined:
        jobs = [(args.combined, samples)]
    else:
        jobs = [(args.prefix + g, s) for g, s in groups.items()]
    # Validate the entire batch before writing, including on case-insensitive filesystems.
    destinations = []
    seen_paths = {"filter_report.csv": "the filter report"}
    for name, cols in jobs:
        filename = safe_name(name) + ".csv"
        key = unicodedata.normalize("NFC", filename).casefold()
        if key in seen_paths:
            ap.error(f"output filename collision: {name!r} and {seen_paths[key]!r} both use {filename!r}")
        seen_paths[key] = name
        destinations.append((os.path.join(args.out_dir, filename), cols))
    os.makedirs(args.out_dir, exist_ok=True)
    for path, cols in destinations:
        block = kept[cols]
        if args.layout == "features":
            write_input(path, list(block.index), block.to_numpy())
        else:
            write_input(path, cols, block.to_numpy().T)
        rows, ncol = (n_kept, len(cols)) if args.layout == "features" else (len(cols), n_kept)
        print(f"wrote {path}: {rows} rows x {ncol} values")
    rep_path = os.path.join(args.out_dir, "filter_report.csv")
    report.to_csv(rep_path, index_label="feature", encoding="utf-8-sig")
    print(f"wrote {rep_path}")


if __name__ == "__main__":
    main()
