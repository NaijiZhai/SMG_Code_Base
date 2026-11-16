# SMG speed test: v1.0 (Excel) vs Python scripts

## Files

- `scripts/real_bench.py`: all steps; see "Run". It times the repository's `cluster.py` and `DAA.py` and reads the maize transcriptome and US adult inputs from the repository's `paper_data/`.
- `excel/SMG v1.0.xlsm`: the original v1.0 workbook.
- `data/pbmc68k_10x/`: 10x Genomics 68k PBMC single-cell data; download it, see "Data".
- `results/`:
  - `real_python_runs.tsv`: Python run times.
  - `real_excel_runs.tsv`: Excel run times and checks (timed on an Excel-saved copy of the same inputs, which differs only in rounding).
  - `real_edges.tsv`: edges below the threshold at each size.
  - `real_extrapolation.tsv`: extrapolated Excel run times and speed-ups.
- Created when the steps run:
  - `inputs/`, `outputs/`.
  - `excel/SMG_v1.0_bench_real.xlsm`: the benchmark workbook; after the Excel run it holds the run times.
  - `excel/SMGBench.bas`: Excel timing module; it answers the macros' input boxes in place of a user and records each run time.
  - `results/real_excel_plan.tsv`: Excel run plan and check values.

## Data

- Maize root transcriptome: Wang D et al. Microbiome 12 (2024): 124.
- US adult gut microbiome: Yatsunenko T et al. Nature 486 (2012): 222–227.
- Both are in the repository's `paper_data/`.
- 10x Genomics 68k PBMC: Zheng GXY et al. Nature Communications 8 (2017): 14049. Download it into `data/pbmc68k_10x/` and leave it compressed (SHA-256 `3f35f37ff344bc9b32f97cd003ac986ebb9b5d7f31006c53dff4cb38da267931`):

```bash
curl -L --create-dirs -o data/pbmc68k_10x/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz https://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz
```

## Run

On macOS, from this folder:

```bash
python -m pip install numpy==2.2.6 pandas==2.2.3 numba==0.61.2 openpyxl==3.1.5
python scripts/real_bench.py prepare    # builds inputs/
python scripts/real_bench.py excel      # builds the benchmark workbook and results/real_excel_plan.tsv
python scripts/real_bench.py python     # times the Python scripts; appends to results/real_python_runs.tsv, so move the old table away before re-timing
```

Excel timing:

1. Open `excel/SMG_v1.0_bench_real.xlsm` in Excel and enable macros.
2. Tools → Macro → Visual Basic Editor, then File → Import File… and import `excel/SMGBench.bas`.
3. Tools → Macro → Macros…, run RunSMGBench. Leave Excel alone until "Finished…" appears.

Summary:

```bash
python scripts/real_bench.py collect        # reads the Excel run times from the workbook
python scripts/real_bench.py edges
python scripts/real_bench.py extrapolate    # needs only the two run-time tables in results/
```
