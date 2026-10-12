#!/usr/bin/env python
"""Step 1d: denoise a harmonized donor screen by projecting each perturbation's delta onto the
leading principal components of all its deltas (gene programmes shared by many perturbations).

A delta from ~50-150 cells is mostly sampling noise (harness split-half cosine ~0.04 on K562 GWPS).
Real knockdown responses share a limited number of programmes; independent noise does not. So:
  1. PCA over the covered genes of all perturbations' deltas, rows weighted by n_cells (noise
     variance ~ 1/n_cells), after removing the weighted mean delta (the screen's shared response);
  2. k = number of eigenvalues above the largest eigenvalue of column-shuffled copies of the
     same matrix (Horn's parallel analysis: shuffling keeps each gene's spread, destroys
     structure), capped at --max-k; or a fixed --n-components;
  3. every perturbation's delta -> mean + projection onto the k components (optionally each
     component shrunk by its signal fraction 1 - noise_edge/eigenvalue, --wiener);
  4. each perturbation's own target gene keeps its measured value (the knockdown is not a shared
     programme; --no-keep-target to drop it). Uncovered genes stay 0.
Output has the same keys as 02_harmonize.py (plus denoise_*), so 09_make_cells.py takes it as a
donor as is.

Usage:
    python scripts/16_denoise_donor.py --npz results/harmonized/replogle_k562_gwps.npz \
        --out results/denoised/replogle_k562_gwps_pca.npz [--n-components K] [--wiener] \
        [--min-cells 50] [--panel data_dir/pert_counts.csv]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse.linalg import eigsh


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True, help="harmonized donor from 02_harmonize.py")
    ap.add_argument("--out", required=True, help="output .npz path")
    ap.add_argument("--n-components", type=int, default=None,
                    help="fixed number of components (default: parallel analysis)")
    ap.add_argument("--max-k", type=int, default=300, help="cap on the parallel-analysis k")
    ap.add_argument("--n-shuffles", type=int, default=3)
    ap.add_argument("--wiener", action="store_true",
                    help="shrink each kept component by 1 - noise_edge/eigenvalue")
    ap.add_argument("--min-cells", type=int, default=0,
                    help="perturbations used to FIT the components (all are projected)")
    ap.add_argument("--no-keep-target", action="store_true")
    ap.add_argument("--panel", default=None,
                    help="CSV with target_gene: also report how much of these perts' deltas survives")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    z = np.load(args.npz, allow_pickle=True)
    data = {k: z[k] for k in z.files}                 # one read per key (each z[k] decompresses)
    perts = data["perts"].astype(str)
    genes = data["genes"].astype(str)
    cov = data["covered"].astype(bool)
    delta = data["delta"].astype(np.float64)
    n_cells = data["n_cells"].astype(float) if "n_cells" in data else np.ones(len(perts))
    cov_idx = np.where(cov)[0]
    D = delta[:, cov_idx]
    fit = n_cells >= args.min_cells
    print(f"{len(perts)} perturbations x {len(cov_idx)} covered genes; fitting on {int(fit.sum())} "
          f"(>= {args.min_cells} cells)", flush=True)

    w = n_cells[fit] / n_cells[fit].mean()
    mu = np.average(D[fit], axis=0, weights=w)
    Xf = (D[fit] - mu) * np.sqrt(w)[:, None]
    nf = Xf.shape[0]
    eig, V = np.linalg.eigh(Xf.T @ Xf / nf)
    order = np.argsort(eig)[::-1]
    eig, V = np.clip(eig[order], 0, None), V[:, order]
    total = eig.sum()

    edges = []
    for s in range(args.n_shuffles):
        Xs = np.take_along_axis(Xf, rng.random(Xf.shape).argsort(axis=0), axis=0)
        edges.append(float(eigsh(Xs.T @ Xs / nf, k=1, return_eigenvectors=False)[0]))
        del Xs
    edge = max(edges)
    k_pa = int((eig > edge).sum())
    k = args.n_components if args.n_components else max(1, min(k_pa, args.max_k))
    print(f"noise edge (max shuffled eigenvalue over {args.n_shuffles}) {edge:.4g}; "
          f"{k_pa} components above it; using k = {k}", flush=True)
    print(f"variance explained by k: {eig[:k].sum() / total:.3f}; "
          f"top eigenvalues {np.round(eig[:5], 4).tolist()}", flush=True)

    Vk = V[:, :k]
    gain = np.clip(1.0 - edge / np.maximum(eig[:k], 1e-30), 0.0, 1.0) if args.wiener else np.ones(k)
    R = D - mu
    Rh = ((R @ Vk) * gain) @ Vk.T
    out = np.zeros_like(delta, dtype=np.float32)
    out[:, cov_idx] = (mu + Rh).astype(np.float32)

    gpos = {g: j for j, g in enumerate(cov_idx)}       # Challenge column -> position in D
    gidx = {g: i for i, g in enumerate(genes)}
    tgt = np.array([gpos.get(gidx.get(p, -1), -1) for p in perts])
    has_t = tgt >= 0
    if not args.no_keep_target:
        rows = np.where(has_t)[0]
        out[rows, cov_idx[tgt[rows]]] = delta[rows, cov_idx[tgt[rows]]]

    # report on deviations from the mean, own target excluded
    Rm, Rhm = R.copy(), Rh.copy()
    rows = np.where(has_t)[0]
    Rm[rows, tgt[rows]] = 0.0
    Rhm[rows, tgt[rows]] = 0.0
    kept = np.linalg.norm(Rhm, axis=1) / np.maximum(np.linalg.norm(Rm, axis=1), 1e-30)
    print(f"fraction of each delta's norm kept (own target excluded): median {np.median(kept):.3f}, "
          f"IQR {np.percentile(kept, 25):.3f}-{np.percentile(kept, 75):.3f}", flush=True)
    by_n = pd.Series(kept).groupby(pd.qcut(n_cells, 3, duplicates="drop")).median().round(3)
    print(f"  by n_cells tertile: {by_n.to_dict()}", flush=True)
    if args.panel:
        panel = set(pd.read_csv(args.panel)["target_gene"].astype(str))
        m = np.isin(perts, list(panel))
        print(f"  panel perts ({int(m.sum())}/{len(panel)} in donor): median kept "
              f"{np.median(kept[m]):.3f}", flush=True)

    data["delta"] = out
    data.update(denoise_k=k, denoise_k_parallel=k_pa, denoise_noise_edge=edge,
                denoise_eigs=eig[: max(k, 50)].astype(np.float32), denoise_wiener=args.wiener,
                denoise_keep_target=not args.no_keep_target, denoise_source=str(args.npz))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, **data)
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
