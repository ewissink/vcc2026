#!/usr/bin/env python
"""Step 1b: reduce one public Perturb-seq h5ad to per-perturbation pseudobulk deltas
on the Challenge gene order (gene_names.csv).

For each perturbation: delta = mean log1p(CPM) of perturbed cells - mean log1p(CPM) of
control cells. Genes absent from the source are zero-filled and flagged in `covered`.

Usage:
    python scripts/02_harmonize.py --h5ad raw/replogle_k562_essential.h5ad \
        --gene-names /path/to/vcc_data/gene_names.csv --name replogle_k562_essential \
        --cell-line K562 --pert-col gene --control-label non-targeting \
        [--gene-col gene_name] [--min-cells 20] --out harmonized/
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def log1p_cpm(X):
    X = sp.csr_matrix(X, dtype=np.float64)
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib == 0] = 1
    Xn = sp.diags(1e6 / lib) @ X
    Xn.data = np.log1p(Xn.data)
    return Xn.tocsr()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--gene-names", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--cell-line", required=True)
    p.add_argument("--pert-col", required=True)
    p.add_argument("--control-label", required=True)
    p.add_argument("--gene-col", default=None, help="var column with symbols; default var_names")
    p.add_argument("--min-cells", type=int, default=20)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    genes = pd.read_csv(args.gene_names, header=None).iloc[:, 0].astype(str)
    if genes.iloc[0].lower() in ("gene", "gene_name", "gene_names", "symbol"):
        genes = genes.iloc[1:]
    genes = genes.reset_index(drop=True)
    gidx = {g: i for i, g in enumerate(genes)}

    a = ad.read_h5ad(args.h5ad)
    src = a.var[args.gene_col].astype(str).values if args.gene_col else a.var_names.astype(str).values
    # keep first occurrence of duplicated symbols
    src_to_target, seen = [], set()
    for j, g in enumerate(src):
        if g in gidx and g not in seen:
            seen.add(g)
            src_to_target.append((j, gidx[g]))
    src_cols = np.array([s for s, _ in src_to_target])
    tgt_cols = np.array([t for _, t in src_to_target])
    covered = np.zeros(len(genes), dtype=bool)
    covered[tgt_cols] = True
    print(f"{len(tgt_cols)}/{len(genes)} challenge genes present in source")

    Xn = log1p_cpm(a.X)[:, src_cols]
    pert = a.obs[args.pert_col].astype(str).values
    is_ctrl = pert == args.control_label
    if is_ctrl.sum() == 0:
        raise SystemExit(f"no cells with {args.pert_col} == {args.control_label!r}")
    ctrl_mean = np.zeros(len(genes), dtype=np.float32)
    ctrl_mean[tgt_cols] = np.asarray(Xn[is_ctrl].mean(axis=0)).ravel()

    df = pd.Series(np.arange(a.n_obs)).groupby(pert).apply(list)
    names, deltas, ncells = [], [], []
    for g, idx in df.items():
        if g == args.control_label or len(idx) < args.min_cells:
            continue
        d = np.zeros(len(genes), dtype=np.float32)
        d[tgt_cols] = np.asarray(Xn[idx].mean(axis=0)).ravel() - ctrl_mean[tgt_cols]
        names.append(g)
        deltas.append(d)
        ncells.append(len(idx))

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out / f"{args.name}.npz",
        perts=np.array(names), delta=np.vstack(deltas), n_cells=np.array(ncells),
        ctrl_mean=ctrl_mean, covered=covered, genes=genes.values,
        cell_line=args.cell_line, n_ctrl=int(is_ctrl.sum()),
    )
    print(f"{args.name}: {len(names)} perturbations kept (>= {args.min_cells} cells)")


if __name__ == "__main__":
    main()
