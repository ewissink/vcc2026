#!/usr/bin/env python
"""Step 1g: build a donor .npz from an HKU DataHub effect_<screen>_moments.h5ad.

The Liu/Wang 2026 release (HKU DataHub 34002315) provides uniformly processed per-perturbation
effect matrices for 12+ CRISPRi screens. Each effect_<ds>_moments.h5ad is:
  rows   = perturbations, obs_names = target-gene symbol (controls are NOT rows)
  cols   = genes, var indexed by gene symbol
  X      = effect E(p,g) = perturbed - control on the log1p(CP10k) scale
  layers = effect_se (and others); uns = control_label, control_profile, batch_column
X is a log1p-normalized mean difference, the same family as 02_harmonize.py's delta (which uses
CP1e6); 09_make_cells divides every donor by ln2, so these average cleanly with our other donors
(no log2 rescaling, unlike the CD4 DESeq2 object). We map to the Challenge gene axis and write the
02_harmonize keys, plus effect_se on the Challenge axis for optional reliability weighting.

Usage:
    python scripts/20_hku_donor.py --h5ad perturb_datasets/hku/effect_Pan-GW-hESC_moments.h5ad \
        --gene-names data_dir/gene_names.csv --cell-line H1 --out results/harmonized/pan_h1.npz
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_genes(path):
    genes = [l.strip().split(",")[0] for l in open(path) if l.strip()]
    if genes and genes[0].lower() in ("gene", "gene_name", "gene_names", "symbol"):
        genes = genes[1:]
    return genes


def dense(M):
    return np.asarray(M.todense()) if sp.issparse(M) else np.asarray(M)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--gene-names", required=True)
    ap.add_argument("--cell-line", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--control-label", default=None, help="override uns['control_label']")
    args = ap.parse_args()

    challenge = read_genes(args.gene_names)
    gidx = {g: i for i, g in enumerate(challenge)}
    G = len(challenge)

    a = ad.read_h5ad(args.h5ad)
    sym = a.var_names.astype(str).values
    src_cols, tgt_cols, seen = [], [], set()
    for j, g in enumerate(sym):
        if g in gidx and g not in seen:
            seen.add(g); src_cols.append(j); tgt_cols.append(gidx[g])
    src_cols, tgt_cols = np.array(src_cols), np.array(tgt_cols)
    covered = np.zeros(G, dtype=bool); covered[tgt_cols] = True
    print(f"{len(tgt_cols)}/{G} challenge genes present in {Path(args.h5ad).name}", flush=True)

    ctrl_label = args.control_label or (str(a.uns["control_label"]) if "control_label" in a.uns else None)
    perts_all = a.obs_names.astype(str).values
    keep_rows = np.array([p != ctrl_label for p in perts_all])      # drop a control row if present
    perts = perts_all[keep_rows]

    X = dense(a.X)[keep_rows]
    X = np.nan_to_num(X.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    delta = np.zeros((len(perts), G), dtype=np.float32)
    delta[:, tgt_cols] = X[:, src_cols].astype(np.float32)

    se_full = np.zeros((len(perts), G), dtype=np.float32)
    if "effect_se" in a.layers:
        S = np.nan_to_num(dense(a.layers["effect_se"])[keep_rows].astype(np.float64), nan=np.inf)
        se_full[:, tgt_cols] = S[:, src_cols].astype(np.float32)

    ctrl_mean = np.zeros(G, dtype=np.float32)
    if "control_profile" in a.uns:
        cp = np.asarray(a.uns["control_profile"], dtype=np.float64).ravel()
        if len(cp) == len(sym):
            ctrl_mean[tgt_cols] = cp[src_cols].astype(np.float32)
    n_cells = (a.obs["n_cells"].astype(float).values[keep_rows]
               if "n_cells" in a.obs.columns else np.zeros(len(perts)))

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        perts=perts.astype(str), delta=delta, covered=covered,
        genes=np.array(challenge, dtype=str), ctrl_mean=ctrl_mean,
        n_cells=np.nan_to_num(n_cells, nan=0.0).astype(float),
        cell_line=args.cell_line, n_ctrl=0,
        effect_se=se_full,
    )
    print(f"{len(perts)} perturbations -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
