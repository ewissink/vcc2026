#!/usr/bin/env python
"""Step 2: PROXY scoring + baselines on the leave-one-cell-line-out folds.

Works on the pseudobulk deltas from 02_harmonize.py (mean log1p-CPM, perturbed - control),
NOT on single cells, so these are rough stand-ins for the vcc2026 metrics, useful for
ranking methods against each other and against baselines. Final checks must use
cell-eval2 on simulated cells.

Proxy metrics, per random 300-perturbation panel from a fold's test perts:
  pds_cosine    close to the official definition: 1 - k_p/(n-1), k_p = (tie-averaged) rank of
                the correct measured effect among all panel effects by cosine distance;
                all panel target genes are removed from the vectors. No-skill = 0.5.
  expr_ratio    sum_p mean_g (pred-meas)^2 / sum_p mean_g meas^2 (no jackknife correction;
                lower is better; predicting no change gives exactly 1).
  dir_fidelity  "significant" = top --top-k genes by |measured delta|. k_match / top_k, where
                k_match = predicted top-k genes that are in the reference set with the same sign.
  dir_reach     reference genes ranked by |pred|; largest prefix with sign purity >= 0.9, / top_k.
  sig_jaccard   Jaccard of predicted top-k and reference top-k gene sets.
  lfc_nmae      sum |pred - meas| / sum |meas| over the reference top-k genes (lower is better).
Each perturbation's own target gene is excluded from every metric. Only genes `covered` in
the held-out dataset are used. Metrics are averaged over perturbations, then over panels.

Methods:
  no_change       predict 0
  train_mean      mean delta over the training perturbations
  oracle_mean     mean measured delta of the panel (the official 'mean-response' baseline b;
                  uses held-out truth, so it is a reference point, not a submittable method)
  transfer        the same gene's delta from the training line(s), n_cells-weighted;
                  falls back to train_mean when the gene was not perturbed in training

Usage:
    python scripts/04_score_baselines.py --harmonized harmonized/ --splits splits.json \
        --fold K562 [--stratum all|novel|seen] [--panel-size 300] [--n-panels 20] [--top-k 100]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def load_datasets(d):
    out = {}
    for f in sorted(Path(d).glob("*.npz")):
        # allow_pickle: older 02_harmonize outputs stored `genes` as an object array.
        # These are our own files, not downloads.
        z = np.load(f, allow_pickle=True)
        out[f.stem] = dict(
            perts=z["perts"].astype(str), delta=z["delta"], n=z["n_cells"].astype(float),
            covered=z["covered"], genes=z["genes"].astype(str), cell_line=str(z["cell_line"]),
        )
    return out


def merge(dsets, names):
    """n_cells-weighted mean delta per perturbation across datasets; returns perts, delta, covered."""
    acc, wsum, cov = {}, {}, None
    for k in names:
        d = dsets[k]
        cov = d["covered"].copy() if cov is None else (cov | d["covered"])
        for i, g in enumerate(d["perts"]):
            acc[g] = acc.get(g, 0) + d["delta"][i] * d["n"][i]
            wsum[g] = wsum.get(g, 0) + d["n"][i]
    perts = np.array(sorted(acc))
    delta = np.vstack([acc[g] / wsum[g] for g in perts]).astype(np.float32)
    return perts, delta, cov


def cosine_rank_scores(P, M):
    """P, M: (n, g). Returns per-row pds score; zero-norm rows give all-tied distances."""
    def unit(A):
        nrm = np.linalg.norm(A, axis=1, keepdims=True)
        return np.where(nrm > 0, A / np.where(nrm > 0, nrm, 1), 0.0)
    D = 1.0 - unit(P) @ unit(M).T                  # distance of pred_i to measured_j
    n = D.shape[0]
    scores = np.empty(n)
    for i in range(n):
        d = D[i]
        dii = d[i]
        tol = 1e-9
        n_less = np.sum(d < dii - tol)
        n_tie = np.sum(np.abs(d - dii) <= tol)
        k = n_less + (n_tie - 1) / 2.0
        scores[i] = 1.0 - k / (n - 1)
    return scores


def score_panel(pred, meas, tgt_idx, top_k):
    """pred, meas: (n, g) over eligible genes; tgt_idx[i]: column of pert i's own target (or -1)."""
    n, g = meas.shape
    panel_targets = [t for t in tgt_idx if t >= 0]
    keep = np.ones(g, dtype=bool)
    keep[panel_targets] = False                    # pds: drop all panel target genes
    pds = cosine_rank_scores(pred[:, keep], meas[:, keep]).mean()

    num = den = 0.0
    fid, reach, jac, nmae = [], [], [], []
    for i in range(n):
        m = meas[i].copy()
        p = pred[i].copy()
        if tgt_idx[i] >= 0:
            m[tgt_idx[i]] = 0.0
            p[tgt_idx[i]] = 0.0
        elig = np.ones(g, dtype=bool)
        if tgt_idx[i] >= 0:
            elig[tgt_idx[i]] = False
        num += np.mean((p[elig] - m[elig]) ** 2)
        den += np.mean(m[elig] ** 2)
        ref = np.argsort(-np.abs(m))[:top_k]
        prd = np.argsort(-np.abs(p))[:top_k]
        match = np.sum(np.isin(prd, ref) & (np.sign(p[prd]) == np.sign(m[prd])))
        fid.append(match / top_k)
        order = ref[np.argsort(-np.abs(p[ref]))]   # reference genes by predicted confidence
        ok = np.cumsum(np.sign(p[order]) == np.sign(m[order])) / np.arange(1, len(order) + 1)
        good = np.where(ok >= 0.9)[0]
        reach.append((good.max() + 1) / top_k if len(good) else 0.0)
        jac.append(len(np.intersect1d(ref, prd)) / len(np.union1d(ref, prd)))
        nmae.append(np.sum(np.abs(p[ref] - m[ref])) / max(np.sum(np.abs(m[ref])), 1e-12))
    return dict(pds_cosine=pds, expr_ratio=num / den, dir_fidelity=np.mean(fid),
                dir_reach=np.mean(reach), sig_jaccard=np.mean(jac), lfc_nmae=np.mean(nmae))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--harmonized", required=True)
    ap.add_argument("--splits", required=True)
    ap.add_argument("--fold", required=True, help="held-out cell line key in splits.json")
    ap.add_argument("--stratum", choices=["all", "novel", "seen"], default="all")
    ap.add_argument("--panel-size", type=int, default=300)
    ap.add_argument("--n-panels", type=int, default=20)
    ap.add_argument("--top-k", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None, help="optional CSV of per-panel results")
    args = ap.parse_args()

    fold = json.loads(Path(args.splits).read_text())[args.fold]
    dsets = load_datasets(args.harmonized)
    genes = next(iter(dsets.values()))["genes"]
    gidx = {g: i for i, g in enumerate(genes)}

    h_perts, h_delta, h_cov = merge(dsets, fold["held_out_datasets"])
    t_perts, t_delta, _ = merge(dsets, fold["train_datasets"])
    t_index = {g: i for i, g in enumerate(t_perts)}
    train_mean = t_delta.mean(axis=0)

    pool = np.array(fold["test_perts"])
    in_train = np.array([g in t_index for g in pool])
    if args.stratum == "novel":
        pool = pool[~in_train]
    elif args.stratum == "seen":
        pool = pool[in_train]
    pool = np.array([g for g in pool if g in set(h_perts)])
    size = min(args.panel_size, len(pool))
    print(f"fold {args.fold}, stratum {args.stratum}: {len(pool)} candidate perts, "
          f"panel size {size}, {int(h_cov.sum())} covered genes")
    if size < 20:
        raise SystemExit("too few perturbations for a meaningful panel")

    cols = np.where(h_cov)[0]
    h_index = {g: i for i, g in enumerate(h_perts)}
    rng = np.random.default_rng(args.seed)
    rows = []
    for panel_i in range(args.n_panels):
        panel = rng.choice(pool, size=size, replace=False)
        meas = np.vstack([h_delta[h_index[g]] for g in panel])[:, cols]
        tgt = np.array([np.where(cols == gidx[g])[0][0] if g in gidx and gidx[g] in set(cols) else -1
                        for g in panel])
        preds = {
            "no_change": np.zeros_like(meas),
            "train_mean": np.tile(train_mean[cols], (size, 1)),
            "oracle_mean": np.tile(meas.mean(axis=0), (size, 1)),
            "transfer": np.vstack([t_delta[t_index[g]] if g in t_index else train_mean
                                   for g in panel])[:, cols],
        }
        for name, P in preds.items():
            rows.append(dict(panel=panel_i, method=name, **score_panel(P, meas, tgt, args.top_k)))

    df = pd.DataFrame(rows)
    if args.out:
        df.to_csv(args.out, index=False)
    summ = df.groupby("method").agg(["mean", "std"]).drop(columns="panel")
    pd.set_option("display.width", 200)
    print(summ.round(3).to_string())


if __name__ == "__main__":
    main()
