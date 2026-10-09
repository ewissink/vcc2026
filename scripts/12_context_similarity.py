#!/usr/bin/env python
"""Rough similarity of each Challenge context (A/B/C) to each donor line, from control profiles.

Compares mean log1p-CPM per gene: 00_inspect_contexts.py's context_means.npz against the
`ctrl_mean` stored in each harmonized donor .npz (02_harmonize.py). Genes used: those covered by
the donor and expressed (mean > --min-mean) in at least one of the two profiles.
Caveat: the two files normalize over different gene sets (all 18,533 genes vs the donor's
covered genes), so absolute correlations are only roughly comparable across donors. Treat the
ranking as a hint, not a measurement.

Usage:
    python scripts/12_context_similarity.py --context-means results/00_inspect/context_means.npz \
        --harmonized results/harmonized [--min-mean 0.1]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--context-means", required=True)
    ap.add_argument("--harmonized", required=True)
    ap.add_argument("--min-mean", type=float, default=0.1)
    args = ap.parse_args()

    cm = np.load(args.context_means, allow_pickle=True)
    genes = cm["genes"].astype(str)
    contexts = [k for k in cm.files if k != "genes"]
    rows = []
    for f in sorted(Path(args.harmonized).glob("*.npz")):
        z = np.load(f, allow_pickle=True)
        if not (z["genes"].astype(str) == genes).all():
            raise SystemExit(f"{f.name}: gene order differs from context_means.npz")
        covered, dmean = z["covered"], z["ctrl_mean"].astype(float)
        for c in contexts:
            cmean = cm[c].astype(float)
            use = covered & ((dmean > args.min_mean) | (cmean > args.min_mean))
            r = np.corrcoef(cmean[use], dmean[use])[0, 1]
            rho = spearmanr(cmean[use], dmean[use])[0]
            rows.append(dict(donor=f.stem, cell_line=str(z["cell_line"]), context=c,
                             n_genes=int(use.sum()), pearson=r, spearman=rho))
    df = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    print(df.round(3).to_string(index=False))
    print("\nSpearman, donors x contexts:")
    print(df.pivot(index="donor", columns="context", values="spearman").round(3))


if __name__ == "__main__":
    main()
