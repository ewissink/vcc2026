#!/usr/bin/env python
"""Step 6: train a cross-cell-line transfer model (ridge) and predict perturbation effects.

Rule-compliant alternative to the donor-blend in 09_make_cells.py: instead of ADDING a scaled
measured donor delta into the prediction, we TRAIN a model on public screens and let the model
predict. The measured deltas are the model's training TARGETS and input FEATURES, never the entry.

Model (per gene g, per (perturbation p, target cell line L)):
  predict  delta_{p,g,L}  from
    * the measured delta of (p, g) in each OTHER donor line  (the transfer signal)
    * that gene's baseline expression (ctrl_mean) in L and in each donor
    * donor coverage flags
  with ridge regression. Weights are SHARED across genes and perturbations (each (p, g) is one
  training row); the features carry all the per-gene / per-line information. Donor columns are
  ordered by cell line, and a perturbation absent from a donor contributes 0 with its flag off,
  so the fitted donor weights are comparable across folds.

Evaluation is leave-one-cell-line-out over splits.json (03_make_splits.py): for each held-out
line we train on all rows from the other lines and predict that line's perturbations. This is the
honest test of cross-line generalization and the number to quote. The single-donor donor-blend of
09 is ~the special case of one donor column with a hand-set weight, so compare against 04.

For a real submission, pass --predict-contexts with a context->donor-line map to emit a predicted
effect table (same gene axis and format 09_make_cells.py consumes via --pred-npz; TODO wire 09).

Each fold also reports a single-donor baseline (each donor's own delta vs the held-out line's,
cosine being scale-free) -- the number the fitted multi-donor model must beat. Features are
standardized before the ridge, and --alphas sweeps the penalty reusing one pass over the data.

Usage (leave-one-line-out proxy evaluation):
    python scripts/18_train_transfer.py --harmonized results/harmonized --splits results/splits.json \
        --alphas 1,10,100 --min-expr 0.5 --out results/transfer_loco.csv
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def load_datasets(d):
    out = {}
    for f in sorted(Path(d).glob("*.npz")):
        z = np.load(f, allow_pickle=True)
        out[f.stem] = dict(
            perts=z["perts"].astype(str), delta=z["delta"].astype(np.float32),
            covered=z["covered"].astype(bool), genes=z["genes"].astype(str),
            ctrl_mean=z["ctrl_mean"].astype(np.float32), cell_line=str(z["cell_line"]),
            n=z["n_cells"].astype(float),
        )
    return out


def line_to_datasets(dsets):
    m = {}
    for k, v in dsets.items():
        m.setdefault(v["cell_line"], []).append(k)
    return m


def merge_line(dsets, names):
    """n_cells-weighted mean delta + ctrl_mean over a line's datasets. Returns {pert: delta}, cov,
    ctrl_mean. All share one gene axis (the Challenge order), checked by the caller."""
    acc, wsum, cov = {}, {}, None
    cm_acc = cm_w = 0.0
    for k in names:
        d = dsets[k]
        cov = d["covered"].copy() if cov is None else (cov | d["covered"])
        w = d["n"].sum()
        cm_acc = cm_acc + d["ctrl_mean"] * w
        cm_w = cm_w + w
        for i, g in enumerate(d["perts"]):
            acc[g] = acc.get(g, 0) + d["delta"][i] * d["n"][i]
            wsum[g] = wsum.get(g, 0) + d["n"][i]
    delta = {g: (acc[g] / wsum[g]).astype(np.float32) for g in acc}
    return delta, cov, (cm_acc / max(cm_w, 1)).astype(np.float32)


def build_rows(target_line, lines, line_delta, line_cov, line_cm, perts, G):
    """Feature matrix X and target y over (perturbation, gene) rows for `perts` in `target_line`.
    donor lines = every line except target_line, in a fixed order."""
    donors = [l for l in lines if l != target_line]
    n = len(perts)
    # per-donor: delta (n, G), coverage (G,), presence of the pert (n,)
    dblocks, cblocks, pblocks, cmblocks = [], [], [], []
    for l in donors:
        d = line_delta[l]
        D = np.zeros((n, G), np.float32)
        present = np.zeros(n, np.float32)
        for i, p in enumerate(perts):
            if p in d:
                D[i] = d[p]
                present[i] = 1.0
        dblocks.append(D)
        cblocks.append(np.broadcast_to(line_cov[l], (n, G)))
        pblocks.append(present)
        cmblocks.append(np.broadcast_to(line_cm[l], (n, G)))
    # features per (p, g): [donor deltas..., donor cov..., target ctrl_mean, donor ctrl_means...,
    #                       donor present...]
    feats = []
    feats += [b.reshape(-1, 1) for b in dblocks]                       # donor delta
    feats += [b.reshape(-1, 1).astype(np.float32) for b in cblocks]    # donor covered
    feats.append(np.broadcast_to(line_cm[target_line], (n, G)).reshape(-1, 1))   # target expr
    feats += [b.reshape(-1, 1) for b in cmblocks]                      # donor expr
    feats += [np.repeat(b, G).reshape(-1, 1) for b in pblocks]         # donor present
    X = np.hstack(feats)
    names = ([f"delta[{l}]" for l in donors] + [f"cov[{l}]" for l in donors]
             + ["expr[target]"] + [f"expr[{l}]" for l in donors] + [f"present[{l}]" for l in donors])
    return X, donors, names


def ridge_precompute(X, y):
    """Center X and y, standardize columns to unit std, and return the small Gram matrices so a
    sweep over alpha reuses one pass over the (large) X. Weights are returned in RAW feature units."""
    mu = X.mean(0)
    sd = X.std(0)
    sd[sd < 1e-8] = 1.0
    Xc = (X - mu) / sd
    ym = float(y.mean())
    XtX = Xc.T @ Xc
    Xty = Xc.T @ (y - ym)
    return dict(mu=mu, sd=sd, ym=ym, XtX=XtX, Xty=Xty, F=Xc.shape[1])


def ridge_solve(pc, alpha):
    """Ridge weights for one alpha from precomputed Gram matrices. Intercept is unpenalized."""
    w_std = np.linalg.solve(pc["XtX"] + alpha * np.eye(pc["F"]), pc["Xty"])
    w = w_std / pc["sd"]                                  # back to raw units
    b = pc["ym"] - pc["mu"] @ w
    return w.astype(np.float32), np.float32(b)


def median_cos(pred, meas, test_perts, elig, gidx):
    """Median cosine(pred_i, meas_i) over eligible genes, each pert's own target dropped."""
    cos = []
    for i, p in enumerate(test_perts):
        e = elig.copy()
        if p in gidx:
            e[gidx[p]] = False
        a, m = pred[i, e], meas[i, e]
        na, nm = np.linalg.norm(a), np.linalg.norm(m)
        if na > 0 and nm > 0:
            cos.append(float(a @ m / (na * nm)))
    return float(np.median(cos)) if cos else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--harmonized", required=True)
    ap.add_argument("--splits", required=True)
    ap.add_argument("--alphas", default="1,10,100", help="comma-separated ridge penalties to sweep")
    ap.add_argument("--min-expr", type=float, default=0.0,
                    help="train/eval only on genes with target-line ctrl_mean above this")
    ap.add_argument("--out", default=None, help="CSV: per-fold weights and proxy cosine")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    dsets = load_datasets(args.harmonized)
    folds = json.load(open(args.splits))
    lines = sorted({v["cell_line"] for v in dsets.values()})
    # one shared gene axis
    genes = next(iter(dsets.values()))["genes"]
    for v in dsets.values():
        assert list(v["genes"]) == list(genes), "all harmonized files must share the gene order"
    G = len(genes)
    gidx = {g: i for i, g in enumerate(genes)}
    l2d = line_to_datasets(dsets)
    line_delta, line_cov, line_cm = {}, {}, {}
    for l in lines:
        line_delta[l], line_cov[l], line_cm[l] = merge_line(dsets, l2d[l])

    alphas = [float(a) for a in args.alphas.split(",")]
    rows = []
    for line in lines:
        if line not in folds:
            continue
        train_lines = [l for l in lines if l != line]
        # TRAIN: stack rows from every (held-out = each training line in turn is the "target"),
        # so the model sees many (donor set -> line) mappings, never the evaluated line.
        Xtr, ytr = [], []
        for tl in train_lines:
            perts_tl = sorted(line_delta[tl])
            X, donors, fnames = build_rows(tl, lines, line_delta, line_cov, line_cm, perts_tl, G)
            y = np.vstack([line_delta[tl][p] for p in perts_tl]).reshape(-1)
            keep = np.broadcast_to(line_cm[tl] > args.min_expr, (len(perts_tl), G)).reshape(-1)
            Xtr.append(X[keep]); ytr.append(y[keep])
        Xtr = np.vstack(Xtr); ytr = np.concatenate(ytr)
        pc = ridge_precompute(Xtr, ytr)
        del Xtr, ytr

        # EVAL on the held-out line's test perts (proxy: cosine of predicted vs measured delta,
        # own target removed), restricted to covered+expressed genes
        test_perts = [p for p in folds[line]["test_perts"] if p in line_delta[line]]
        X, donors, fnames = build_rows(line, lines, line_delta, line_cov, line_cm, test_perts, G)
        meas = np.vstack([line_delta[line][p] for p in test_perts])
        elig = (line_cm[line] > args.min_expr) & line_cov[line]

        # single-donor baseline: cosine is scale-free, so this is just each donor's own delta vs
        # the held-out line's delta (no fitting). The number multi-donor must beat.
        di = {name: j for j, name in enumerate(donors)}
        base = {}
        for d in donors:
            dcol = X[:, di[d]].reshape(len(test_perts), G)    # that donor's delta block
            base[d] = median_cos(dcol, meas, test_perts, elig, gidx)
        best_d = max(base, key=base.get)
        print(f"{line}: test perts {len(test_perts)}; best single donor "
              f"{best_d} cos {base[best_d]:.3f} (all: "
              f"{ {k: round(v,3) for k,v in base.items()} })", flush=True)

        for alpha in alphas:
            w, b = ridge_solve(pc, alpha)
            pred = (X @ w + b).reshape(len(test_perts), G)
            mc = median_cos(pred, meas, test_perts, elig, gidx)
            wd = dict(zip(fnames, np.round(w, 4)))
            dw = {k: v for k, v in wd.items() if k.startswith('delta[')}
            print(f"   alpha {alpha:>6g}: multi-donor cos {mc:.3f}  "
                  f"vs best-single {base[best_d]:.3f}   weights {dw}", flush=True)
            rows.append(dict(line=line, n_test=len(test_perts), alpha=alpha, median_cos=mc,
                             best_single_donor=best_d, best_single_cos=base[best_d],
                             intercept=float(b), **wd))

    if args.out:
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
