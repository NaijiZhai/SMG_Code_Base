import csv
import os
import statistics
import subprocess
import sys
import tarfile
import time
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GENES = ROOT.parent / 'paper_data/Corn Root Microbiome and Transcriptome/Transcriptome'
REF = ROOT.parent / 'paper_data/Human Gut Microbiome/Result 3/3DAA/Adult/input_3DAA_U_Adult.csv'
PBMC = ROOT / 'data/pbmc68k_10x/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz'
V1 = ROOT / 'excel/SMG v1.0.xlsm'
IN, OUT, RES, XL = ROOT / 'inputs', ROOT / 'outputs', ROOT / 'results', ROOT / 'excel'

# python runs: the first k genes (24 samples) and the first n cells (536 genes)
G_SIZES = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 16249]
C_SIZES = [6 * 2 ** k for k in range(14)] + [68579]
# excel runs: only the sizes expected to finish within about 15 minutes
EXCEL_RUNS = ([('Graph', 'genes', k) for k in (10, 20, 40, 80, 160, 320, 640)]
              + [('Graph', 'cells', n) for n in (6, 12, 24)]
              + [('SymDiff', 'genes', k) for k in (10, 20, 40, 80, 160, 320)])


def read_table(path):
    # headerless csv: ID, then values
    rows = list(csv.reader(open(path, encoding='utf-8-sig')))
    return [r[0] for r in rows], np.array([[float(v) for v in r[1:]] for r in rows])


def write_table(path, ids, values):
    # headerless csv, 0 for zero, other values at full precision
    with open(path, 'w') as f:
        for name, row in zip(ids, values):
            f.write(name + ',' + ','.join('0' if v == 0 else repr(float(v)) for v in row) + '\n')


def read_daa(path):
    # DAA output: total changed edges, then each OTU's count
    rows = [r for r in csv.reader(open(path, encoding='utf-8-sig')) if r]
    return int(rows[1][1]), {r[0]: int(r[1]) for r in rows[2:]}


def load_cluster():
    # cluster.py imported for its functions
    # numba's cache from an imported copy breaks the command-line runs of the same file: keep it elsewhere
    os.environ['NUMBA_CACHE_DIR'] = str(OUT / 'numba_cache_import')
    spec = importlib.util.spec_from_file_location('cluster', ROOT.parent / 'cluster.py')
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def tables():
    # sorted tables written by prepare
    gw, gl = np.load(IN / 'genes_wt.npy'), np.load(IN / 'genes_lateral.npy')
    cells = np.load(IN / 'cells.npy', mmap_mode='r')
    return (gw, gl, (IN / 'genes_ids.txt').read_text().split()), (cells, (IN / 'cells_ids.txt').read_text().split())


def prepare():
    # output folders
    for d in (IN, OUT, RES, XL):
        d.mkdir(exist_ok=True)

    # genes: WT and lateral-root mutants, genes ranked by the larger of the two group means
    ids, A = read_table(GENES / 'input_DAA_GWT.csv')
    ids2, B = read_table(GENES / 'input_DAA_G34.csv')
    assert ids == ids2
    order = np.argsort(-np.maximum(A.mean(axis=1), B.mean(axis=1)), kind='stable')
    np.save(IN / 'genes_wt.npy', A[order])
    np.save(IN / 'genes_lateral.npy', B[order])
    (IN / 'genes_ids.txt').write_text('\n'.join(ids[i] for i in order) + '\n')
    print('genes:', A.shape, B.shape, 'zero share', round(np.mean(A == 0), 4), round(np.mean(B == 0), 4))

    # cells: 10x 68k PBMC read from the archive, matrix.mtx starts with % comment lines
    with tarfile.open(PBMC, 'r:gz') as tar:
        f = tar.extractfile('filtered_matrices_mex/hg19/matrix.mtx')
        line = f.readline()
        while line.startswith(b'%'):
            line = f.readline()
        n_genes, n_cells, _ = (int(v) for v in line.split())
        t = pd.read_csv(f, sep=' ', header=None, names=['g', 'c', 'v'])
        names = [r.split('\t')[0] for r in tar.extractfile('filtered_matrices_mex/hg19/genes.tsv').read().decode().splitlines()]
    # each cell's UMIs / its total, the 536 genes with the highest mean, cells shuffled once
    g, c, v = t.g.to_numpy() - 1, t.c.to_numpy() - 1, t.v.to_numpy(float)
    rel = v / np.bincount(c, weights=v, minlength=n_cells)[c]
    top = np.argsort(-np.bincount(g, weights=rel, minlength=n_genes), kind='stable')[:536]
    rank = np.full(n_genes, -1)
    rank[top] = np.arange(536)
    keep = rank[g] >= 0
    M = np.zeros((536, n_cells))
    M[rank[g[keep]], c[keep]] = rel[keep]
    M = M[:, np.random.default_rng(61115).permutation(n_cells)]
    np.save(IN / 'cells.npy', M)
    (IN / 'cells_ids.txt').write_text('\n'.join(names[i] for i in top) + '\n')
    print('cells:', M.shape, 'zero share', round(np.mean(M == 0), 4))

    # reference: US adults, the 536 OTUs with the highest mean
    ids, X = read_table(REF)
    top = np.argsort(-X.mean(axis=1), kind='stable')[:536]
    write_table(IN / 'ref_us_adult.csv', [ids[i] for i in top], X[top])
    print('ref:', (536, X.shape[1]), 'zero share', round(np.mean(X[top] == 0), 4))


VBA = r'''Attribute VB_Name = "SMGBench"
Option Explicit
' runs the unchanged macros Graph and SymDiff on the sizes in sheet Bench_plan and times them
' InputBox and MsgBox below replace VBA's own in this workbook: during a run they answer the prompts from the plan

Private active As Boolean
Private answers As Variant
Private nextAnswer As Long

Public Function InputBox(Prompt As Variant, Optional Title As Variant, Optional Default As Variant, _
        Optional XPos As Variant, Optional YPos As Variant, Optional HelpFile As Variant, Optional Context As Variant) As String
    If Not active Then
        InputBox = VBA.InputBox(Prompt)
    Else
        InputBox = CStr(answers(nextAnswer))
        nextAnswer = nextAnswer + 1
    End If
End Function

Public Function MsgBox(Prompt As Variant, Optional Buttons As Variant, Optional Title As Variant, _
        Optional HelpFile As Variant, Optional Context As Variant) As VbMsgBoxResult
    If Not active Then MsgBox = VBA.MsgBox(Prompt) Else MsgBox = vbOK
End Function

Public Sub RunSMGBench()
    Dim wb As Workbook, plan As Worksheet, res As Worksheet
    Dim r As Long, m As Long, n As Long, n2 As Long, kind As String, t0 As Double, sec As Double
    Dim mat As Variant, i As Long, j As Long, total As Double, le As Long, p1 As Double
    Set wb = ThisWorkbook
    Set plan = wb.Worksheets("Bench_plan")
    Set res = wb.Worksheets("Bench_results")
    res.Cells.Clear
    res.Range("A1:F1").Value = Array("run_id", "seconds", "matrix1_sum", "edges1_le", "compare_edges", "excel")
    For r = 2 To plan.Cells(plan.Rows.Count, 1).End(xlUp).Row
        kind = plan.Cells(r, 2).Value
        m = plan.Cells(r, 5).Value
        n = plan.Cells(r, 6).Value
        wb.Worksheets("Original1").Cells.Clear
        wb.Worksheets("Original2").Cells.Clear
        wb.Worksheets("Original1").Range("A1").Resize(m, n + 1).Value = wb.Worksheets(plan.Cells(r, 3).Value).Range("A1").Resize(m, n + 1).Value
        If kind = "SymDiff" Then
            n2 = plan.Cells(r, 9).Value
            wb.Worksheets("Original2").Range("A1").Resize(m, n2 + 1).Value = wb.Worksheets(plan.Cells(r, 4).Value).Range("A1").Resize(m, n2 + 1).Value
            answers = Array("1", "0.1", CStr(plan.Cells(r, 7).Value), CStr(plan.Cells(r, 8).Value))
        Else
            answers = Array(CStr(plan.Cells(r, 7).Value))
        End If
        wb.Worksheets("Start").Activate
        nextAnswer = 0
        active = True
        t0 = Timer
        If kind = "Graph" Then Module1.Graph Else Module2.SymDiff
        sec = Timer - t0
        If sec < 0 Then sec = sec + 86400
        active = False
        ' checks, not timed: distance sum and edges with distance <= threshold; SymDiff's changed edges
        mat = wb.Worksheets("Matrix1").Range("B2").Resize(m, m).Value
        p1 = CDbl(plan.Cells(r, 7).Value)
        total = 0: le = 0
        For i = 1 To m
            For j = 1 To m
                total = total + mat(i, j)
                If j > i And mat(i, j) <= p1 Then le = le + 1
            Next j
        Next i
        res.Cells(r, 1).Value = plan.Cells(r, 1).Value
        res.Cells(r, 2).Value = sec
        res.Cells(r, 3).Value = total
        res.Cells(r, 4).Value = le
        If kind = "SymDiff" Then res.Cells(r, 5).Value = wb.Worksheets("Compare").Cells(1, 2).Value
        wb.Save
    Next r
    res.Range("F2").Value = Application.Version & " | " & Application.OperatingSystem
    wb.Save
    VBA.MsgBox "Finished. The timings are in the sheet Bench_results."
End Sub
'''


def excel():
    # run plan: python's threshold and check values for each excel run
    import openpyxl
    cl = load_cluster()
    (gw, gl, gids), (cells, cids) = tables()
    data = {'genes': ('Bench_g1', gids, gw, 'Bench_g2', gl), 'cells': ('Bench_c', cids, cells, '', None)}
    plan = []
    for run_id, (kind, series, size) in enumerate(EXCEL_RUNS, 1):
        sheet1, ids, A, sheet2, B = data[series]
        A = A[:size] if series == 'genes' else np.asarray(A[:, :size])
        a = cl.prepare_values(A)
        t1 = cl.group_threshold(a, 90.0, 1.0, None)
        d1 = cl.distance_matrix(a)
        iu = np.triu_indices(len(a), 1)
        row = {'run_id': run_id, 'kind': kind, 'data1': sheet1, 'data2': '', 'otus': len(A), 'samples': A.shape[1],
               'p1': format(t1, '.17g'), 'p2': '', 'samples2': '', 'matrix1_sum': d1.sum(), 'edges1_le': int((d1[iu] <= t1).sum()),
               'changed_lt': ''}
        if kind == 'SymDiff':
            b = cl.prepare_values(B[:size])
            t2 = cl.group_threshold(b, 90.0, 1.0, None)
            d2 = cl.distance_matrix(b)
            # v1.0 links pairs with distance < bound
            row.update(data2=sheet2, p2=format(t2, '.17g'), samples2=b.shape[1],
                       changed_lt=int(((d1[iu] < t1) != (d2[iu] < t2)).sum()))
        plan.append(row)
    plan = pd.DataFrame(plan)
    plan.to_csv(RES / 'excel_plan.tsv', sep='\t', index=False)

    # benchmark workbook: the original plus the input sheets, the run plan and an empty results sheet
    wb = openpyxl.load_workbook(V1, keep_vba=True)
    for title, ids, values in (('Bench_g1', gids[:640], gw[:640]), ('Bench_g2', gids[:640], gl[:640]),
                               ('Bench_c', cids, np.asarray(cells[:, :24]))):
        ws = wb.create_sheet(title)
        for name, row in zip(ids, values):
            ws.append([name, *[float(v) for v in row]])
    ws = wb.create_sheet('Bench_plan')
    cols = ['run_id', 'kind', 'data1', 'data2', 'otus', 'samples', 'p1', 'p2', 'samples2']
    ws.append(cols)
    for r in plan[cols].itertuples(index=False):
        ws.append(list(r))
    # thresholds as text, so all 17 digits reach the macros
    for row in ws.iter_rows(min_row=2, min_col=7, max_col=8):
        for cell in row:
            cell.data_type = 's'
    wb.create_sheet('Bench_results')
    wb.save(XL / 'SMG_v1.0_bench.xlsm')
    # timing module to import into the workbook, crlf line endings
    (XL / 'SMGBench.bas').write_bytes(VBA.replace('\n', '\r\n').encode('ascii'))
    print(plan[['run_id', 'kind', 'otus', 'samples', 'samples2', 'p1', 'edges1_le', 'changed_lt']].to_string(index=False))


def python():
    # time cluster.py and DAA.py from the command line: start-up, reading, computing and writing
    py = sys.executable
    (gw, gl, gids), (cells, cids) = tables()
    out = RES / 'python_runs.tsv'

    def run(cmd):
        # wall time and peak memory of one run, ru_maxrss in bytes on macOS
        t0 = time.perf_counter()
        p = subprocess.Popen([str(x) for x in cmd], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        _, status, usage = os.wait4(p.pid, 0)
        assert status == 0, cmd
        return time.perf_counter() - t0, usage.ru_maxrss

    def timed(series, tool, paths, rows, samples, samples2=''):
        # repeats set by the first run: 5 under 2 s, 3 under 20 s, otherwise 1; the median is kept
        cmd = ([py, ROOT.parent / 'cluster.py', '-i', paths[0], '-o', OUT / 'distances.csv'] if tool == 'cluster' else
               [py, ROOT.parent / 'DAA.py', '-in1', paths[0], '-in2', paths[1], '-o', OUT / 'daa.csv'])
        wall, peak = run(cmd)
        walls = [wall]
        while len(walls) < (5 if walls[0] < 2 else 3 if walls[0] < 20 else 1):
            wall, rss = run(cmd)
            walls.append(wall)
            peak = max(peak, rss)
        row = {'series': series, 'tool': tool, 'rows': rows, 'samples': samples, 'samples2': samples2, 'runs': len(walls),
               'seconds_median': round(statistics.median(walls), 4), 'seconds_min': round(min(walls), 4),
               'peak_rss_gb': round(peak / 1e9, 3), 'daa_edges': '', 'same_as_package': ''}
        if tool == 'daa':
            row['daa_edges'], deg = read_daa(OUT / 'daa.csv')
            # all genes: the manuscript's comparison, checked against output_DAA_GWT_G34.csv
            if rows == len(gids):
                row['same_as_package'] = (row['daa_edges'], deg) == read_daa(GENES / 'output_DAA_GWT_G34.csv')
        # one row per size, appended to the table
        new = not out.exists()
        with open(out, 'a', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(row), delimiter='\t')
            if new:
                w.writeheader()
            w.writerow(row)
        print('\t'.join(str(v) for v in row.values()), flush=True)

    def gene_inputs(k):
        # the first k genes of both groups as input files
        paths = [IN / f'g{k}_wt.csv', IN / f'g{k}_lateral.csv']
        write_table(paths[0], gids[:k], gw[:k])
        write_table(paths[1], gids[:k], gl[:k])
        return paths

    # first runs compile numba's functions (cached afterwards): not timed
    p = gene_inputs(10)
    run([py, ROOT.parent / 'cluster.py', '-i', p[0], '-o', OUT / 'distances.csv'])
    run([py, ROOT.parent / 'DAA.py', '-in1', p[0], '-in2', p[1], '-o', OUT / 'daa.csv'])
    # genes: cluster.py on WT, DAA.py on WT vs lateral-root mutants
    for k in G_SIZES:
        paths = gene_inputs(k)
        timed('genes', 'cluster', paths, k, gw.shape[1])
        timed('genes', 'daa', paths, k, gw.shape[1], gl.shape[1])
    # cells: cluster.py on the first n cells
    for n in C_SIZES:
        path = IN / f'c{n}.csv'
        write_table(path, cids, np.asarray(cells[:, :n]))
        timed('cells', 'cluster', [path], 536, n)
    # reference: cluster.py on US adults
    timed('ref', 'cluster', [IN / 'ref_us_adult.csv'], 536, 136)


def edges():
    # edges (pairs with distance < threshold) at each size: Graph cuts and inserts a row and a column for each one,
    # so they set the reordering part of its run time; excel holds at most 16,383 samples
    cl = load_cluster()
    (gw, gl, gids), (cells, cids) = tables()
    jobs = ([('genes', gw[:k]) for k in G_SIZES] + [('cells', np.asarray(cells[:, :n])) for n in C_SIZES if n <= 16383]
            + [('ref', read_table(IN / 'ref_us_adult.csv')[1])])
    rows = []
    for series, X in jobs:
        a = cl.prepare_values(X)
        t = cl.group_threshold(a, 90.0, 1.0, None)
        d = cl.distance_matrix(a)
        m = len(a)
        # symmetric matrix with a zero diagonal
        e = int(((d < t).sum() - m) // 2)
        rows.append({'series': series, 'rows': m, 'samples': X.shape[1], 'threshold': t, 'pairs': m * (m - 1) // 2,
                     'edges_lt': e, 'edge_share': e / (m * (m - 1) / 2)})
        print('\t'.join(str(v) for v in rows[-1].values()), flush=True)
        del d
    pd.DataFrame(rows).to_csv(RES / 'edges.tsv', sep='\t', index=False)


def collect():
    # excel run times and checks from the workbook, next to python's check values
    import openpyxl
    rows = list(openpyxl.load_workbook(XL / 'SMG_v1.0_bench.xlsm', read_only=True, data_only=True)['Bench_results'].values)
    res = pd.DataFrame([r[:5] for r in rows[1:] if r[0] is not None], columns=list(rows[0][:5]))
    both = pd.read_csv(RES / 'excel_plan.tsv', sep='\t').merge(res, on='run_id', suffixes=('', '_excel'))
    both['matrix1_rel_diff'] = (both.matrix1_sum_excel - both.matrix1_sum).abs() / both.matrix1_sum
    both['excel'] = rows[1][5]
    both.to_csv(RES / 'excel_runs.tsv', sep='\t', index=False)
    print(both[['run_id', 'kind', 'otus', 'samples', 'seconds', 'matrix1_rel_diff', 'edges1_le', 'edges1_le_excel',
                'changed_lt', 'compare_edges']].to_string(index=False))


def extrapolate():
    # excel has 16,384 columns and v1.0 keeps the names in column A
    LIMIT = 16383
    py = pd.read_csv(RES / 'python_runs.tsv', sep='\t')
    xl = pd.read_csv(RES / 'excel_runs.tsv', sep='\t')
    # v1.0 beyond its largest run: strictly with the square of the size and, for genes, also with the slope of its
    # three largest runs on log-log axes
    rows = []
    for name, kind, data1, tool, col in (('genes Graph', 'Graph', 'Bench_g1', 'cluster', 'otus'),
                                         ('genes SymDiff', 'SymDiff', 'Bench_g1', 'daa', 'otus'),
                                         ('cells Graph', 'Graph', 'Bench_c', 'cluster', 'samples')):
        series = name.split()[0]
        e = xl[(xl.kind == kind) & (xl.data1 == data1)].sort_values(col)
        size, sec = e[col].to_numpy(float), e.seconds.to_numpy()
        slope = np.polyfit(np.log(size[-3:]), np.log(sec[-3:]), 1)[0] if series == 'genes' else np.nan
        # reference point (US adults, 536 x 136): from the largest cell run, square only
        if series == 'cells':
            reference_seconds = sec[-1] * (136 / size[-1]) ** 2
        p = py[(py.series == series) & (py.tool == tool)]
        for s, t in zip(p.rows if series == 'genes' else p.samples, p.seconds_median):
            if s > LIMIT:
                continue
            r = {'comparison': name, 'rows': s if series == 'genes' else 536,
                 'samples': ('24 vs 23' if kind == 'SymDiff' else 24) if series == 'genes' else s,
                 'growth_exponent': slope, 'python_seconds': t,
                 'excel_measured': sec[size == s][0] if s in size else np.nan,
                 'excel_square': sec[-1] * (s / size[-1]) ** 2 if s > size[-1] else np.nan,
                 'excel_growth': sec[-1] * (s / size[-1]) ** slope if s > size[-1] else np.nan}
            rows.append(r)
    ref = py[py.series == 'ref'].iloc[0]
    rows.append({'comparison': 'ref Graph', 'rows': 536, 'samples': 136, 'growth_exponent': np.nan,
                 'python_seconds': ref.seconds_median, 'excel_measured': np.nan,
                 'excel_square': reference_seconds, 'excel_growth': np.nan})
    out = pd.DataFrame(rows)
    # speed-ups: excel time / python time
    for c in ('measured', 'square', 'growth'):
        out['speedup_' + c] = out['excel_' + c] / out.python_seconds
    out.to_csv(RES / 'extrapolation.tsv', sep='\t', index=False)
    print(out.to_string(index=False))


if __name__ == '__main__':
    {'prepare': prepare, 'excel': excel, 'python': python, 'edges': edges, 'collect': collect,
     'extrapolate': extrapolate}[sys.argv[1]]()
