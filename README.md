# Standalone SMG scripts

Four Python scripts for the Slope-Matrix-Graph (SMG) algorithm. Each script
runs on its own: download the one you need, install its dependencies and run it
with your input files.

| Script | Purpose | Output |
|---|---|---|
| `clean.py` | Convert counts to relative abundances and filter OTUs | One input CSV per group and a filter report |
| `cluster.py` | Build one OTU network and find its connected components | OTU distance matrix and cluster member list (CSV) |
| `DAA.py` | Compare two OTU networks | Changed-edge counts (CSV) |
| `BDA.py` | Compute sample-to-sample SMG distances | Sample distance matrix (CSV) |

## Installation

Tested with Python 3.12. Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

The analysis scripts use NumPy, Pandas and Numba. The cleaner uses NumPy and
Pandas, plus openpyxl when reading an Excel input. Run any script with `-h` to
see its options.

## Cleaning

The count table has OTU IDs in its first column and sample IDs in its header.
The group table has columns named `sample` and `group`.

```bash
python clean.py -c counts.csv -g groups.csv -o cleaned
```

The cleaner converts sample columns to relative abundances, then keeps an OTU
when its largest group mean relative abundance is at least 0.00001 (0.001%) and
its largest group prevalence is at least 0.2 (20%). It writes headerless
`input_<group>.csv` files and `filter_report.csv`.

Filter all groups used in DAA or BDA together, so that their input files contain
the same OTUs. To cluster one group, use a group table containing only that group.
For BDA, `--layout samples` writes one sample per row and `--combined`
writes all samples into one file:

```bash
python clean.py -c counts.csv -g groups.csv -o cleaned_bda --layout samples --combined all_samples
python BDA.py -i cleaned_bda/all_samples.csv -o output/sample_distances.csv
```

## Analysis input

`cluster.py`, `DAA.py` and `BDA.py` read headerless CSV files: the first column
is the ID, the other columns are abundances. `cluster.py` and `DAA.py` take one
OTU per row, `BDA.py` one sample per row. IDs must be unique and abundances
must be non-negative numbers. Output paths must end in `.csv`.

## Clustering

```bash
python cluster.py -i cleaned/input_A.csv -o output/otu_distances.csv
```

This writes two files:

- `output/otu_distances.csv`: the distance matrix of all OTUs, ordered by
  connected component (largest first), with isolated OTUs last.
- `output/otu_distances_clusters.csv`: the components with at least two OTUs,
  one OTU per row (`cluster_id`, `cluster_size`, `otu_id`), in the same order
  as the matrix. The second row contains the column headers. To read this file with
Pandas, use `pd.read_csv(path, skiprows=1)`.

## Differential network analysis

```bash
python DAA.py -in1 cleaned/input_A.csv -in2 cleaned/input_B.csv -o output/daa.csv
```

The output CSV has a `Direct Comparison` row, an `Edges` row with the number of
edges present in one network but not the other, and then each OTU with its
number of such edges, largest first. OTUs are matched between the two files by
ID.

## Network threshold

Clustering and DAA min-max scale each OTU separately and, by default, use the
group-adaptive arc threshold at 90 degrees. `-a`/`--angle` sets the arc angle,
`-l`/`--lam` multiplies the threshold, and `--column-normalize` rescales each
sample to relative abundance before the min-max scaling. `-t` (clustering) or
`-pos1`/`-pos2` (DAA) set a fixed threshold.

## Beta difference analysis

```bash
python BDA.py -i samples.csv -o output/sample_distances.csv
```

Each row is a sample, with the OTUs in the same column order for every sample.
BDA computes SMG distances on these abundances as given and writes the sample
distance matrix in the input sample order.

## Reproducing the manuscript

`paper_data/` holds the data behind the manuscript's results. In each dataset
folder, `Data.xlsx` has the raw counts, one sheet per group, and `groups/` has
one group table per analysis, with the groups that are cleaned together. The
steps below rebuild all other files in place in about a minute; run them from
the repository folder. Intermediate files go to `output/`.

Put the sheets of each dataset side by side in one count table:

```bash
mkdir -p output
python - <<'EOF'
import csv
import openpyxl

def merge(book, sheets, out):
    # every sheet lists the same features in the same order; the baby sheets start with a row of ages
    parts = []
    for ws in openpyxl.load_workbook(book, read_only=True, data_only=True).worksheets:
        if ws.title.startswith(sheets):
            rows = [r for r in ws.iter_rows(values_only=True) if any(v is not None for v in r)]
            rows = rows[1:] if rows[0][0] == 'Age' else rows
            cols = [j for j, v in enumerate(rows[0]) if j and v is not None]
            parts.append([[str(r[0])] + [r[j] for j in cols] for r in rows])
    assert len({len(p) for p in parts}) == 1
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        for i, rows in enumerate(zip(*parts)):
            assert i == 0 or len({r[0] for r in rows}) == 1
            w.writerow([rows[0][0]] + [v for r in rows for v in r[1:]])

merge('paper_data/Human Gut Microbiome/Data.xlsx', '', 'output/gut_counts.csv')
merge('paper_data/Corn Root Microbiome and Transcriptome/Data.xlsx', ('Rhizosphere', 'Root'), 'output/otu_counts.csv')
merge('paper_data/Corn Root Microbiome and Transcriptome/Data.xlsx', 'Transcriptome', 'output/gene_counts.csv')
EOF
```

Clean the groups of each analysis together and put the input files in place:

```bash
G="paper_data/Human Gut Microbiome"
M="paper_data/Corn Root Microbiome and Transcriptome"
for n in Amerindians Malawians US; do
  c=${n:0:1}
  python clean.py -c output/gut_counts.csv -g "$G/groups/Result1_$n.csv" -o output/1BDA_$c --layout samples --combined input_1BDA_$c
  # AmzC1babyF1.418864 has no age: it is cleaned with the babies but left out of the file
  grep -v '^AmzC1babyF1\.418864,' output/1BDA_$c/input_1BDA_$c.csv > "$G/Result 1/$n/input_1BDA_$c.csv"
  python clean.py -c output/gut_counts.csv -g "$G/groups/Result2_$n.csv" -o output/2BDA_$c --layout samples --prefix input_2BDA_
  cp output/2BDA_$c/input_*.csv "$G/Result 2/$n/"
done
for a in Baby Child Adult; do
  python clean.py -c output/gut_counts.csv -g "$G/groups/Result3_$a.csv" -o output/3BDA_$a --layout samples --combined input_3BDA_$a
  cp output/3BDA_$a/input_3BDA_$a.csv "$G/Result 3/3BDA/$a/"
  python clean.py -c output/gut_counts.csv -g "$G/groups/Result3_$a.csv" -o output/3DAA_$a --prefix input_3DAA_
  cp output/3DAA_$a/input_*.csv "$G/Result 3/3DAA/$a/"
done
python clean.py -c output/otu_counts.csv -g "$M/groups/Rhizosphere.csv" -o output/Rhizosphere --prefix input_DAA_
python clean.py -c output/otu_counts.csv -g "$M/groups/Root.csv" -o output/Root --prefix input_DAA_
python clean.py -c output/gene_counts.csv -g "$M/groups/Transcriptome.csv" -o output/Transcriptome --prefix input_DAA_
for t in Rhizosphere Root Transcriptome; do cp output/$t/input_*.csv "$M/$t/"; done
```

Run BDA and DAA:

```bash
G="paper_data/Human Gut Microbiome"
M="paper_data/Corn Root Microbiome and Transcriptome"
for f in "$G/Result 1"/*/input_*.csv "$G/Result 2"/*/input_*.csv "$G/Result 3/3BDA"/*/input_*.csv; do
  python BDA.py -i "$f" -o "${f/input_/output_}"
done
for a in Baby Child Adult; do
  d="$G/Result 3/3DAA/$a"
  python DAA.py -in1 "$d/input_3DAA_A_$a.csv" -in2 "$d/input_3DAA_M_$a.csv" -o "$d/output_3DAA_A_M_$a.csv"
  python DAA.py -in1 "$d/input_3DAA_A_$a.csv" -in2 "$d/input_3DAA_U_$a.csv" -o "$d/output_3DAA_A_U_$a.csv"
  python DAA.py -in1 "$d/input_3DAA_M_$a.csv" -in2 "$d/input_3DAA_U_$a.csv" -o "$d/output_3DAA_M_U_$a.csv"
done
d="$M/Rhizosphere"
python DAA.py -in1 "$d/input_DAA_SWT.csv" -in2 "$d/input_DAA_S34.csv" -o "$d/output_DAA_SWT_S34.csv"
python DAA.py -in1 "$d/input_DAA_SWT.csv" -in2 "$d/input_DAA_S567.csv" -o "$d/output_DAA_SWT_S567.csv"
d="$M/Root"
python DAA.py -in1 "$d/input_DAA_RWT.csv" -in2 "$d/input_DAA_R34.csv" -o "$d/output_DAA_RWT_R34.csv"
python DAA.py -in1 "$d/input_DAA_RWT.csv" -in2 "$d/input_DAA_R567.csv" -o "$d/output_DAA_RWT_R567.csv"
d="$M/Transcriptome"
python DAA.py -in1 "$d/input_DAA_GWT.csv" -in2 "$d/input_DAA_G34.csv" -o "$d/output_DAA_GWT_G34.csv"
python DAA.py -in1 "$d/input_DAA_GWT.csv" -in2 "$d/input_DAA_G567.csv" -o "$d/output_DAA_GWT_G567.csv"
```

Summarize:

```bash
python - <<'EOF'
import csv
import numpy as np
import pandas as pd

G, M = 'paper_data/Human Gut Microbiome', 'paper_data/Corn Root Microbiome and Transcriptome'
C = {'A': 'Amerindians', 'M': 'Malawians', 'U': 'US'}
AGES = ('Baby', 'Child', 'Adult')
PAIRS = (('A', 'M'), ('A', 'U'), ('M', 'U'))

def read(path, encoding='utf-8'):
    with open(path, encoding=encoding) as f:
        return list(csv.reader(f))

def write(path, rows):
    with open(path, 'w', newline='') as f:
        csv.writer(f, lineterminator='\n').writerows(rows)

def within(m):
    # mean over the sample pairs of one group
    d = m.to_numpy()
    return d[np.triu_indices_from(d, 1)].mean()

def edges(path):
    # second row of a DAA output: the number of changed edges
    return int(read(path, 'utf-8-sig')[1][1])

# Figure 2: each baby's mean distance to the adults of its country, then the adults' own mean
rows = [('country', 'sample', 'age', 'average_distance_to_adults')]
for c, n in C.items():
    groups = read(f'{G}/groups/Result1_{n}.csv')[1:]
    m = pd.read_csv(f'{G}/Result 1/{n}/output_1BDA_{c}.csv', index_col=0)
    adults = [s for s, g, _ in groups if g == 'Adult']
    rows += [(n, s, age, m.loc[s, adults].to_numpy().mean()) for s, g, age in groups if g == 'Baby' and s in m.index]
    rows.append((n, 'Adult_In_Group_Average', '', within(m.loc[adults, adults])))
write(f'{G}/Result 1/Result_1BDA.csv', rows)
# Table 1: mean distance within each age group; its percentages are the child and adult means over the baby mean
write(f'{G}/Result 2/Result_2BDA.csv', [('country', *AGES)] + [
    (n, *[within(pd.read_csv(f'{G}/Result 2/{n}/output_2BDA_{c}_{a}.csv', index_col=0)) for a in AGES])
    for c, n in C.items()])
# Tables 2A and 2B: mean distance and changed edges between the countries of one age group
head = ('age', *[f'{C[x]} vs {C[y]}' for x, y in PAIRS])
t2a, t2b = [head], [head]
for a in AGES:
    m = pd.read_csv(f'{G}/Result 3/3BDA/{a}/output_3BDA_{a}.csv', index_col=0)
    s = {c: [x for x, g in read(f'{G}/groups/Result3_{a}.csv')[1:] if g == f'{c}_{a}'] for c in C}
    t2a.append((a, *[m.loc[s[x], s[y]].to_numpy().mean() for x, y in PAIRS]))
    t2b.append((a, *[edges(f'{G}/Result 3/3DAA/{a}/output_3DAA_{x}_{y}_{a}.csv') for x, y in PAIRS]))
write(f'{G}/Result 3/3BDA/Result_3BDA.csv', t2a)
write(f'{G}/Result 3/3DAA/Result_3DAA.csv', t2b)
# Table 4: changed edges between B73 and each mutant group
tissues = (('Rhizosphere', 'S'), ('Root', 'R'), ('Transcriptome', 'G'))
write(f'{M}/Result_DAA.csv', [('comparison', *[d for d, _ in tissues])] + [
    (label, *[edges(f'{M}/{d}/output_DAA_{t}WT_{t}{g}.csv') for d, t in tissues])
    for label, g in (('B73 vs rum1 & lrt1', '34'), ('B73 vs rth3 & rth5 & rth6', '567'))])
EOF
```

| Summary | Manuscript |
|---|---|
| `Human Gut Microbiome/Result 1/Result_1BDA.csv` | Figure 2 |
| `Human Gut Microbiome/Result 2/Result_2BDA.csv` | Table 1 |
| `Human Gut Microbiome/Result 3/3BDA/Result_3BDA.csv` | Table 2A |
| `Human Gut Microbiome/Result 3/3DAA/Result_3DAA.csv` | Table 2B |
| `Corn Root Microbiome and Transcriptome/Result_DAA.csv` | Table 4 |


## Sources

The original SMG algorithm:
Zhang M, Li X, Oladeinde A, Rothrock M Jr, Pokoo-Aikins A, Zock G.
A Novel Slope-Matrix-Graph Algorithm to Analyze Compositional Microbiome Data.
Microorganisms 12 (2024): 1866. https://doi.org/10.3390/microorganisms12091866

Data:
Wang D et al. Lateral root enriched Massilia associated with plant flowering in
maize. Microbiome 12 (2024): 124. https://doi.org/10.1186/s40168-024-01839-4

Yatsunenko T et al. Human gut microbiome viewed across age and geography.
Nature 486 (2012): 222-227. https://doi.org/10.1038/nature11053
