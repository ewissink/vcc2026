#!/usr/bin/env python
"""Rank DepMap/CCLE cell lines by similarity to each Challenge context (and to known donors).

Query profiles:
  * contexts A/B/C: mean log1p-CPM per gene from 00_inspect_contexts.py (context_means.npz)
  * donors (optional): `ctrl_mean` from harmonized .npz files, as a sanity check -- if the method
    cannot put K562/HepG2/Jurkat near the top of the CCLE ranking, trust the context results less.
Reference: a DepMap expression table (rows = cell lines, columns = genes named "SYMBOL (ID)",
e.g. OmicsExpressionProteinCodingGenesTPMLogp1.csv) plus optionally Model.csv for names/lineages.
Spearman correlation over the most variable genes across CCLE lines (rank-based, so the
different units -- single-cell log1p-CPM vs bulk log2(TPM+1) -- do not matter).
CCLE is bulk RNA-seq and the queries are single-cell, so expect modest correlations and
near-ties between related lines; read the ranking as a short list, not a single answer.
File/column names are from memory of recent DepMap releases and not checked against a download.

Usage:
    python scripts/13_identify_context.py --expression OmicsExpressionProteinCodingGenesTPMLogp1.csv \
        --model-csv Model.csv --context-means results/00_inspect/context_means.npz \
        --harmonized results/harmonized [--n-var-genes 5000] [--top 8]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

EXPECT = {"K562": "K562", "HepG2": "HEPG2", "Jurkat": "JURKAT", "RPE1": "RPE1"}


def unit_ranks(M):
    R = np.apply_along_axis(rankdata, -1, M)
    R = R - R.mean(axis=-1, keepdims=True)
    return R / np.linalg.norm(R, axis=-1, keepdims=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expression", required=True)
    ap.add_argument("--model-csv", default=None)
    ap.add_argument("--context-means", required=True)
    ap.add_argument("--harmonized", default=None)
    ap.add_argument("--n-var-genes", type=int, default=5000)
    ap.add_argument("--top", type=int, default=8)
    args = ap.parse_args()

    ref = pd.read_csv(args.expression, index_col=0)
    ref.columns = [c.split(" (")[0] for c in ref.columns]
    ref = ref.loc[:, ~pd.Index(ref.columns).duplicated()]
    names = pd.DataFrame(index=ref.index)
    names["name"] = ref.index.astype(str)
    names["lineage"] = ""
    if args.model_csv:
        m = pd.read_csv(args.model_csv).set_index("ModelID")
        for src, dst in (("StrippedCellLineName", "name"), ("OncotreeLineage", "lineage")):
            if src in m.columns:
                names[dst] = m[src].reindex(ref.index).fillna(names[dst] if dst == "name" else "").values
    print(f"CCLE reference: {ref.shape[0]} lines x {ref.shape[1]} genes", flush=True)

    cm = np.load(args.context_means, allow_pickle=True)
    genes = cm["genes"].astype(str)
    queries = {f"context {k}": (pd.Series(cm[k].astype(float), index=genes), None)
               for k in cm.files if k != "genes"}
    if args.harmonized:
        for f in sorted(Path(args.harmonized).glob("*.npz")):
            z = np.load(f, allow_pickle=True)
            line = str(z["cell_line"])
            if line in EXPECT and f"donor {line}" not in queries:
                s = pd.Series(z["ctrl_mean"].astype(float), index=z["genes"].astype(str))
                queries[f"donor {line}"] = (s[z["covered"].astype(bool)], EXPECT[line])

    sd = ref.std(axis=0)
    for qname, (prof, expect) in queries.items():
        common = [g for g in prof.index if g in ref.columns]
        top_genes = sd[common].sort_values(ascending=False).index[: args.n_var_genes]
        R = unit_ranks(ref[top_genes].to_numpy(float))
        q = unit_ranks(prof[top_genes].to_numpy(float)[None, :])[0]
        rho = R @ q
        order = np.argsort(-rho)
        print(f"\n== {qname}  ({len(top_genes)} genes)")
        for i in order[: args.top]:
            print(f"  {rho[i]:.3f}  {names['name'].iloc[i]:<18} {names['lineage'].iloc[i]}")
        if expect:
            hit = [r for r, i in enumerate(order, 1) if expect in str(names["name"].iloc[i]).upper()]
            print(f"  expected '{expect}': " + (f"rank {hit[0]} of {len(order)}" if hit else "not in reference"))


if __name__ == "__main__":
    main()
