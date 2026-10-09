#!/usr/bin/env python
"""Smoke test: does cell-eval2's vcc2026 preset run on CPU, and does `pred` need control cells?

Builds a toy `real` (controls + perturbed cells) and a toy `pred` (a noisy copy of the real
effects), then scores pred-with-controls and pred-without-controls. Prints each result or
the error. Run on BioHPC in an env with `pip install cell-eval2` (Python >= 3.11).

Usage:  python scripts/10_celleval_smoke.py [--device cpu]
"""
import argparse
import dataclasses
import traceback

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def make_toy(seed=0, n_genes=300, n_perts=20, n_cells=100, n_ctrl=600):
    rng = np.random.default_rng(seed)
    base = rng.uniform(2, 30, n_genes)                       # mean counts, mostly > 5 CPM-ish
    genes = [f"G{i}" for i in range(n_genes)]
    perts = [f"G{i}" for i in range(n_perts)]
    effects = np.zeros((n_perts, n_genes))
    for i in range(n_perts):
        idx = rng.choice(n_genes, 25, replace=False)
        effects[i, idx] = rng.normal(0, 1.0, 25)             # log2 fc

    def cells(mean, n):
        return rng.negative_binomial(5, 5 / (5 + mean), size=(n, len(mean)))

    def build(effect_noise, include_ctrl):
        X, lab = [], []
        if include_ctrl:
            X.append(cells(base, n_ctrl)); lab += ["non-targeting"] * n_ctrl
        for i, p in enumerate(perts):
            e = effects[i] + rng.normal(0, effect_noise, n_genes) * (effect_noise > 0)
            X.append(cells(base * 2.0 ** e, n_cells)); lab += [p] * n_cells
        a = ad.AnnData(sp.csr_matrix(np.vstack(X).astype(np.float32)),
                       obs=pd.DataFrame({"target_gene": lab}))
        a.obs_names = [f"c{i}" for i in range(a.n_obs)]
        a.var_names = genes
        return a

    return build(0.0, True), build(0.3, True), build(0.3, False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()

    from cell_eval2 import EvalConfig, compute_metrics

    real, pred_with, pred_without = make_toy()
    real.write_h5ad("toy_real.h5ad")
    pred_with.write_h5ad("toy_pred_with_ctrl.h5ad")
    pred_without.write_h5ad("toy_pred_no_ctrl.h5ad")
    cfg = dataclasses.replace(EvalConfig.from_preset("vcc2026"), pert_col="target_gene",
                              control="non-targeting", device=args.device)
    for name in ("toy_pred_no_ctrl.h5ad", "toy_pred_with_ctrl.h5ad"):
        print(f"\n=== pred = {name} (device={args.device})", flush=True)
        try:
            print(compute_metrics(name, "toy_real.h5ad", config=cfg))
        except Exception:
            traceback.print_exc()


if __name__ == "__main__":
    main()
