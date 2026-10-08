#!/usr/bin/env python
"""Step 0: QC the three validation control files and guess which cell line each is.

Usage:
    python scripts/00_inspect_contexts.py --data-dir /path/to/vcc_data --out results/00_inspect

Writes:
    context_summary.csv      per-context cell/library-size/sparsity stats
    marker_expression.csv    mean log1p(CPM) of marker genes per context
    context_correlation.csv  Pearson r of mean log1p(CPM) between contexts
    context_means.npz        mean log1p(CPM) per gene per context (reused downstream)
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Rough lineage markers; this is a heuristic aid, not a classifier.
MARKERS = {
    "K562": ["HBG1", "HBG2", "GYPA", "KIT", "GATA1", "ALAS2"],
    "HepG2": ["ALB", "APOA1", "APOB", "SERPINA1", "AHSG", "TTR"],
    "Jurkat": ["CD3D", "CD3E", "LCK", "IL2RG", "CD2"],
    "RPE1": ["TERT", "RPE65", "BEST1", "VIM", "CRYAB"],
    "hESC": ["POU5F1", "NANOG", "SOX2", "LIN28A", "DPPA4"],
    "HEK293": ["NEFM", "NEFL", "PAX2", "LHX2", "GATA3"],
    "A549": ["SFTPB", "AKR1C1", "AKR1B10", "KRT19", "TFPI2"],
    "MCF7": ["ESR1", "TFF1", "GATA3", "PGR", "KRT18"],
    "HCT116": ["CEACAM5", "KRT20", "CDX2", "LGALS4"],
    "HeLa": ["KRT17", "KRT5", "CDKN2A", "S100A2"],
    "epithelial": ["EPCAM", "KRT8", "KRT18", "CDH1"],
    "mesenchymal": ["VIM", "FN1", "CDH2", "ZEB1"],
}


def log1p_cpm_mean(X):
    """Mean over cells of log1p(CPM); X is raw counts (sparse or dense)."""
    X = sp.csr_matrix(X)
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib == 0] = 1
    scale = sp.diags(1e6 / lib)
    Xn = scale @ X
    Xn.data = np.log1p(Xn.data)
    return np.asarray(Xn.mean(axis=0)).ravel()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    data = Path(args.data_dir)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(data / "gene_names.csv", header=None).iloc[:, 0].astype(str)
    if genes.iloc[0].lower() in ("gene", "gene_name", "gene_names", "symbol"):
        genes = genes.iloc[1:].reset_index(drop=True)
    manifest = json.loads((data / "manifest.json").read_text())
    print("manifest:", json.dumps(manifest, indent=2)[:1500])

    perts = pd.read_csv(data / "pert_counts.csv")
    print(f"{len(perts)} perturbations to predict; columns: {list(perts.columns)}")

    summary, means = [], {}
    for path in sorted(data.glob("context_*.h5ad")):
        ctx = path.stem.split("_", 1)[1]
        a = ad.read_h5ad(path)
        X = sp.csr_matrix(a.X)
        lib = np.asarray(X.sum(axis=1)).ravel()
        order_ok = list(a.var_names) == list(genes)
        is_int = bool(np.all(np.mod(X.data, 1) == 0))
        row = dict(
            context=ctx,
            n_cells=a.n_obs,
            n_genes=a.n_vars,
            gene_order_matches=order_ok,
            raw_integer_counts=is_int,
            median_umi=float(np.median(lib)),
            mean_umi=float(lib.mean()),
            frac_zero=1 - X.nnz / (X.shape[0] * X.shape[1]),
            genes_detected_per_cell=float(np.median(np.diff(X.indptr))),
            n_ntc_guides=a.obs["ntc_id"].nunique() if "ntc_id" in a.obs else np.nan,
        )
        summary.append(row)
        means[ctx] = log1p_cpm_mean(X)
        print(row)

    summ = pd.DataFrame(summary)
    summ.to_csv(out / "context_summary.csv", index=False)

    mean_df = pd.DataFrame(means, index=list(a.var_names))
    np.savez(out / "context_means.npz", genes=np.array(mean_df.index), **{k: mean_df[k].values for k in mean_df})

    rows = []
    for line, mk in MARKERS.items():
        for g in mk:
            if g in mean_df.index:
                rows.append(dict(group=line, gene=g, **mean_df.loc[g].to_dict()))
    pd.DataFrame(rows).to_csv(out / "marker_expression.csv", index=False)

    mean_df.corr().to_csv(out / "context_correlation.csv")
    print("\nbetween-context correlation of mean log1p(CPM):")
    print(mean_df.corr().round(3))
    print("\nmarker expression (mean log1p CPM):")
    print(pd.DataFrame(rows).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
