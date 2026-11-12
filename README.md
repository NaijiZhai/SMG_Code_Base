# SMG Algorithm Implementation

Python implementation of the updated **Slope-Matrix-Graph (SMG) algorithm** for compositional microbiome data. The
algorithm was first published in _Microorganisms_ [[1]](https://doi.org/10.3390/microorganisms12091866). 

| Script | What it does |
|---|---|
| `smg_clean.py` | Converts a count table to relative abundances, filters the OTUs and writes one input file per group |
| `DAA-11-05.py` | Differential network analysis (DAA): counts the edges present in one group's OTU network but not in the other's |
| `cluster-11-05.py` | Builds the OTU network of one group and reports its clusters |
| `BDA-09-03.py` | Beta difference analysis (BDA): computes the SMG distance between every two samples |

Keep all `.py` files in one folder. Run a script with `-h` to see all of its options.

## Installation

Tested with Python 3.12.7.

```bash
git clone https://github.com/NaijiZhai/SMG_Code_Base
cd SMG_Code_Base
pip install -r requirement.txt
```

## Usage

```bash
python smg_clean.py -c counts.csv -g groups.csv -o cleaned
python DAA-11-05.py -in1 cleaned/input_A.csv -in2 cleaned/input_B.csv -o daa.csv
python smg_clean.py -c counts.csv -g groups_A.csv -o cleaned_A
python cluster-11-05.py -i cleaned_A/input_A.csv -o clusters_A.xlsx
python BDA-09-03.py -i samples.csv -o distances.xlsx
```

> [!IMPORTANT]
> **For large data sets, `--format binary` is highly recommended.** Clustering and BDA write Excel workbooks by
> default, which hold every pairwise distance and become very slow and memory-hungry with thousands of OTUs or
> samples. Binary output writes the network and clusters (clustering) or the distance matrix (BDA) to a result
> directory in a fraction of the time.

- `smg_clean.py` takes a count table (OTU IDs in the first column, sample IDs in the header row) and a group table
  (columns `sample` and `group`). It writes one `input_<group>.csv` per group, with one row per OTU and no header, and
  a `filter_report.csv`. All groups to be compared are filtered together (`groups.csv`); a group to be clustered is
  filtered on its own (`groups_A.csv` lists only group A). `--layout samples` writes one row per sample, as BDA needs.
- `DAA-11-05.py` writes the number of changed edges, then every OTU with the number of changed edges it belongs to.
- `cluster-11-05.py` writes an Excel workbook with the distance matrix colored by cluster, plus the `Threshold` and
  `Standard OTU` sheets. For large networks, `--format binary -o clusters_A` writes a folder with `edges.csv`,
  `clusters.csv` and threshold diagnostics. Add `--distances` in binary mode to also save `distances.npy` and `ids.csv`.
- `BDA-09-03.py` writes the distance matrix between the samples to an Excel workbook. For large sample sets,
  `--format binary -o bda_result` is faster: it writes `distances.npy`, with the sample order in `ids.csv`.
- Input CSVs for DAA, clustering and BDA have no header row. The first column contains OTU IDs for DAA and clustering,
or sample IDs for BDA. All BDA samples must use the same OTU set in the same column order.

For binary output, use a separate result directory for each analysis.

## Example of use

The `examples/` folder holds the maize OTU table of Wang et al. [[2]](https://doi.org/10.1186/s40168-024-01839-4)
with its group tables (`WT`: B73; `34`: lateral-root mutants; `567`: root-hair mutants), and the gut data of the
Amerindian babies, children and adults of Yatsunenko et al. [[3]](https://doi.org/10.1038/nature11053). 

**Differential network analysis** of the root samples, B73 vs. lateral-root mutants:

```bash
mkdir -p examples/output
python smg_clean.py -c examples/raw/wang_otu_counts.csv -g examples/groups/root.csv -o examples/output/root
python DAA-11-05.py -in1 examples/output/root/input_WT.csv -in2 examples/output/root/input_34.csv -o examples/output/DAA_WT_vs_34.csv
```

`DAA_WT_vs_34.csv` reports `Edges,5202,` on its second line: 5,202 changed edges, the Root value in Table 4.
With `input_567.csv` as the second file the result is 1,991; the rhizosphere samples
(`examples/groups/rhizosphere.csv`) give 1,620 and 701.

**Clustering** of the B73 root samples, filtered on their own:

```bash
python smg_clean.py -c examples/raw/wang_otu_counts.csv -g examples/groups/root_WT.csv -o examples/output/root_WT
python cluster-11-05.py -i examples/output/root_WT/input_WT.csv -o examples/output/cluster_WT.xlsx
```

The 536 OTUs that pass the filter form 21 clusters with at least two OTUs.

**Beta difference analysis** of the 32 Amerindian babies, and their mean distance:

```bash
python BDA-09-03.py -i examples/gut/amerindian_baby.csv -o examples/output/BDA_amerindian_baby.xlsx
python -c "import numpy as np, pandas as pd; d = pd.read_excel('examples/output/BDA_amerindian_baby.xlsx', index_col=0).to_numpy(); print(f'{d[np.triu_indices_from(d, 1)].mean():.3f}')"
```

This prints 4292.229, the value for Amerindian babies in Table 1; the child and adult files give 2,959.760 and 3,305.305.

## Citation

If you use this implementation in your research, please cite the original paper:

> Zhang M, Li X, Oladeinde A, Rothrock M Jr, Pokoo-Aikins A, Zock G. A Novel Slope-Matrix-Graph Algorithm to Analyze
> Compositional Microbiome Data. _Microorganisms_ 12, no. 9 (2024): 1866.
> [https://doi.org/10.3390/microorganisms12091866](https://doi.org/10.3390/microorganisms12091866)

The example data are from:

> Wang D, He X, Baer M, Lami K, Yu B, Tassinari A, Salvi S, Schaaf G, Hochholdinger F, Yu P. Lateral root enriched
> _Massilia_ associated with plant flowering in maize. _Microbiome_ 12 (2024): 124.
> [https://doi.org/10.1186/s40168-024-01839-4](https://doi.org/10.1186/s40168-024-01839-4)

> Yatsunenko T, Rey FE, Manary MJ, et al. Human gut microbiome viewed across age and geography. _Nature_ 486 (2012):
> 222–227. [https://doi.org/10.1038/nature11053](https://doi.org/10.1038/nature11053)
