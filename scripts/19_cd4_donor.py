#!/usr/bin/env python
"""Step 1f: build a donor .npz from the Marson CD4+ T-cell genome-scale Perturb-seq DE object.

The CZI/Marson release (GWCD4i.DE_stats.h5ad) is already per-perturbation DESeq2 stats, not cells:
obs = 33,983 (perturbation x culture_condition) rows, var = 10,282 genes (var_names Ensembl,
symbols in var['gene_name']), and layers['log_fc'] = log2 fold change vs non-targeting. We take
one culture condition (default Rest, closest to the Jurkat/context-A T-cell line), map to the
Challenge gene axis, and write the same keys as 02_harmonize.py so 09/04/17 consume it directly.

Units: 09_make_cells divides a donor's `delta` by ln2 to get log2. The other donors store
(mean log1p-CPM perturbed - control) ~ natural-log units; this object is already log2 FC, so we
store delta = log_fc * ln2 and 09's /ln2 recovers the true log2 FC. Directionally comparable to
the other donors for averaging; magnitudes are a cleaner DESeq2 estimate.

Also carries per-pert reliability hints from the DE obs (ontarget_effect_size, n_cells_target,
guide_correlation_all) for optional 17-style weighting later.

Usage:
    python scripts/19_cd4_donor.py --de-h5ad perturb_datasets/cd4/GWCD4i.DE_stats.h5ad \
        --gene-names data_dir/gene_names.csv --condition Rest --cell-line CD4T \
        --out results/harmonized/cd4_rest.npz
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def read_genes(path):
    genes = [l.strip().split(",")[0] for l in open(path) if l.strip()]
    if genes and genes[0].lower() in ("gene", "gene_name", "gene_names", "symbol"):
        genes = genes[1:]
    return genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de-h5ad", required=True)
    ap.add_argument("--gene-names", required=True)
    ap.add_argument("--condition", default="Rest")
    ap.add_argument("--cell-line", default="CD4T")
    ap.add_argument("--pert-col", default="target_contrast_gene_name")
    ap.add_argument("--cond-col", default="culture_condition")
    ap.add_argument("--layer", default="log_fc", help="layer holding log2 fold change")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    challenge = read_genes(args.gene_names)
    gidx = {g: i for i, g in enumerate(challenge)}
    G = len(challenge)

    a = ad.read_h5ad(args.de_h5ad, backed="r")
    sym = (a.var["gene_name"].astype(str).values if "gene_name" in a.var.columns
           else a.var_names.astype(str).values)
    # source column j -> Challenge index, first occurrence of each symbol
    src_cols, tgt_cols, seen = [], [], set()
    for j, g in enumerate(sym):
        if g in gidx and g not in seen:
            seen.add(g); src_cols.append(j); tgt_cols.append(gidx[g])
    src_cols, tgt_cols = np.array(src_cols), np.array(tgt_cols)
    covered = np.zeros(G, dtype=bool); covered[tgt_cols] = True
    print(f"{len(tgt_cols)}/{G} challenge genes present in CD4 object", flush=True)

    cond = a.obs[args.cond_col].astype(str).values
    rows = np.where(cond == args.condition)[0]
    if len(rows) == 0:
        raise SystemExit(f"no rows with {args.cond_col} == {args.condition!r}; "
                         f"have {sorted(set(cond))}")
    sub = a[rows].to_memory()
    lfc = np.asarray(sub.layers[args.layer], dtype=np.float64)       # (n_rows, 10282), log2 FC
    lfc = np.nan_to_num(lfc, nan=0.0, posinf=0.0, neginf=0.0)
    pert = sub.obs[args.pert_col].astype(str).values

    # one delta per perturbation (average if a gene appears more than once in this condition)
    def obs_num(col):
        return (pd.to_numeric(sub.obs[col], errors="coerce").values
                if col in sub.obs.columns else np.full(len(sub), np.nan))
    eff = obs_num("ontarget_effect_size"); ncell = obs_num("n_cells_target")
    gcorr = obs_num("guide_correlation_all")

    order = pd.unique(pert)
    delta = np.zeros((len(order), G), dtype=np.float32)
    rel_eff = np.zeros(len(order)); rel_nc = np.zeros(len(order)); rel_gc = np.zeros(len(order))
    LN2 = np.log(2.0)
    for i, p in enumerate(order):
        m = np.where(pert == p)[0]
        v = lfc[m].mean(axis=0)                                     # genes x (source order)
        delta[i, tgt_cols] = (v[src_cols] * LN2).astype(np.float32)  # log2 -> our convention
        rel_eff[i] = np.nanmean(eff[m]); rel_nc[i] = np.nanmean(ncell[m]); rel_gc[i] = np.nanmean(gcorr[m])
    print(f"{len(order)} perturbations in condition {args.condition}", flush=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        perts=order.astype(str), delta=delta, covered=covered,
        genes=np.array(challenge, dtype=str),
        ctrl_mean=np.zeros(G, dtype=np.float32),                    # not used as a donor
        n_cells=np.nan_to_num(rel_nc, nan=0.0).astype(float),
        cell_line=args.cell_line, n_ctrl=0,
        cd4_ontarget_effect=rel_eff.astype(np.float32),
        cd4_guide_correlation=rel_gc.astype(np.float32),
        cd4_condition=args.condition,
    )
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
