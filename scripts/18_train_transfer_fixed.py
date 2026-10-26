"""Corrected LOCO diagnostic. Proxy metrics only; no compliance determination."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def load_datasets(d):
    out = {}
    for f in sorted(Path(d).glob("*.npz")):
        z = np.load(f, allow_pickle=True)
        required = {"perts", "delta", "covered", "genes",
                    "ctrl_mean", "cell_line", "n_cells"}
        if "perts" not in z.files:
            print(f"Skipping non-donor archive {f.name}: keys={z.files}",
                  flush=True)
            z.close()
            continue
        missing = required - set(z.files)
        if missing:
            z.close()
            raise ValueError(f"{f}: missing donor fields {sorted(missing)}")
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
    """Merge effects with coverage-aware weights and explicit missing-count handling."""
    acc, wsum = {}, {}
    G = len(dsets[names[0]]["genes"])
    cov = np.zeros(G, dtype=bool)
    cm_acc = np.zeros(G, dtype=np.float64)
    cm_w = np.zeros(G, dtype=np.float64)

    for k in names:
        d = dsets[k]
        covered = np.asarray(d["covered"], dtype=bool)
        counts = np.asarray(d["n"], dtype=float)
        if len(counts) != len(d["perts"]):
            raise ValueError(f"{k}: cell-count and perturbation lengths differ")

        valid_counts = np.isfinite(counts) & (counts > 0)
        weights = np.where(valid_counts, counts, 1.0)
        if not valid_counts.all():
            examples = d["perts"][~valid_counts][:5].tolist()
            print(
                f"{k}: {int((~valid_counts).sum())}/{len(counts)} perts "
                f"have missing/zero/invalid cell counts; using weight 1. "
                f"Examples: {examples}", flush=True)

        cov |= covered
        cm = np.asarray(d["ctrl_mean"], dtype=float)
        cm_ok = covered & np.isfinite(cm)
        screen_weight = float(weights.sum())
        cm_acc[cm_ok] += cm[cm_ok] * screen_weight
        cm_w[cm_ok] += screen_weight

        for i, pert in enumerate(d["perts"]):
            values = np.asarray(d["delta"][i], dtype=np.float64)
            if not np.isfinite(values[covered]).all():
                raise ValueError(
                    f"{k}, perturbation {pert}: non-finite source effects "
                    "on covered genes; inspect this archive before proceeding")
            if pert not in acc:
                acc[pert] = np.zeros(G, dtype=np.float64)
                wsum[pert] = np.zeros(G, dtype=np.float64)
            acc[pert][covered] += values[covered] * weights[i]
            wsum[pert][covered] += weights[i]

    delta = {}
    for pert in acc:
        np.divide(acc[pert], wsum[pert], out=acc[pert],
                  where=wsum[pert] > 0)
        delta[pert] = acc[pert].astype(np.float32)
    cm = np.divide(cm_acc, cm_w, out=np.zeros(G), where=cm_w > 0)
    return delta, cov, cm.astype(np.float32)


def line_similarity(lines, line_cm, min_expr):
    """Pearson r between every pair of lines' control profiles (ctrl_mean), over genes expressed
    (> min_expr) in BOTH. Computable from control cells alone, so it is available for a held-out
    test line from its NTC profile. Returns {(a, b): r}."""
    sim = {}
    for a in lines:
        for b in lines:
            if a == b:
                continue
            m = (line_cm[a] > min_expr) & (line_cm[b] > min_expr)
            if m.sum() < 2:
                sim[(a, b)] = 0.0
            else:
                sim[(a, b)] = float(np.corrcoef(line_cm[a][m], line_cm[b][m])[0, 1])
    return sim


def build_rows(target_line, lines, line_delta, line_cov, line_cm, perts, G, sim=None):
    """Feature matrix X and target y over (perturbation, gene) rows for `perts` in `target_line`.
    donor lines = every line except target_line, in a fixed order. With `sim`, also add each
    donor's delta times similarity(target_line, donor) so the model can weight donors by how much
    the target line resembles them (selection), not just a global constant (averaging)."""
    donors = list(lines)
    n = len(perts)
    dblocks, cblocks, pblocks, cmblocks = [], [], [], []
    for l in donors:
        d = line_delta[l] if l != target_line else {}
        D = np.zeros((n, G), np.float32)
        present = np.zeros(n, np.float32)
        for i, p in enumerate(perts):
            if p in d:
                D[i] = d[p]
                present[i] = 1.0
        dblocks.append(D)
        cblocks.append(np.broadcast_to(line_cov[l] if l != target_line else np.zeros(G, bool), (n, G)))
        pblocks.append(present)
        cmblocks.append(np.broadcast_to(line_cm[l] if l != target_line else np.zeros(G), (n, G)))
    feats = []
    feats += [b.reshape(-1, 1) for b in dblocks]                       # donor delta
    feats += [b.reshape(-1, 1).astype(np.float32) for b in cblocks]    # donor covered
    feats.append(np.broadcast_to(line_cm[target_line], (n, G)).reshape(-1, 1))   # target expr
    feats += [b.reshape(-1, 1) for b in cmblocks]                      # donor expr
    feats += [np.repeat(b, G).reshape(-1, 1) for b in pblocks]         # donor present
    names = ([f"delta[{l}]" for l in donors] + [f"cov[{l}]" for l in donors]
             + ["expr[target]"] + [f"expr[{l}]" for l in donors] + [f"present[{l}]" for l in donors])
    if sim is not None:
        # delta x similarity(target, donor): sim is constant within one fold, so this is
        # identifiable only because sim VARIES across the stacked training folds.
        feats += [(dblocks[j] * np.float32(sim.get((target_line, l), 0.0))).reshape(-1, 1)
                  for j, l in enumerate(donors)]
        names += [f"delta_x_sim[{l}]" for l in donors]
    X = np.hstack(feats)
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
        cos.append(float(a @ m / (na * nm)) if na > 0 and nm > 0 else 0.0)
    return float(np.median(cos)) if cos else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--harmonized", required=True)
    ap.add_argument("--splits", required=True)
    ap.add_argument("--alphas", default="1,10,100", help="comma-separated ridge penalties to sweep")
    ap.add_argument("--min-expr", type=float, default=0.0,
                    help="train/eval only on genes with target-line ctrl_mean above this")
    ap.add_argument("--out", default=None, help="CSV: per-fold weights and proxy cosine")
    ap.add_argument("--no-sim", action="store_true",
                    help="disable similarity-conditioned donor weights (reproduces global averaging)")
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

    sim = None if args.no_sim else line_similarity(lines, line_cm, args.min_expr)
    if sim is not None:
        print("control-profile similarity (target <- donor), per line:")
        for a in lines:
            print(f"  {a}: " + ", ".join(f"{b} {sim[(a,b)]:.2f}" for b in lines if b != a), flush=True)
    alphas = [float(a) for a in args.alphas.split(",")]
    rows = []
    for line in lines:
        if line not in folds:
            continue
        train_lines = [l for l in lines if l != line]
        # TRAIN: stack rows from every (held-out = each training line in turn is the "target"),
        # so the model sees many (donor set -> line) mappings, never the evaluated line.
        # Approximately equal training weight per target line.
        # Accumulate sufficient statistics without retaining the full design matrix.
        sx = sy = sxx = sxy = None
        total = 0.0
        for tl in train_lines:
            perts_tl = sorted(line_delta[tl])
            nr = len(perts_tl) * int((line_cm[tl] > args.min_expr).sum())
            if nr == 0:
                print(f"Skipping training target {tl}: no eligible rows", flush=True)
                continue
            for lo in range(0, len(perts_tl), 32):
                pp = perts_tl[lo:lo + 32]
                X, donors, fnames = build_rows(
                    tl, train_lines, line_delta, line_cov, line_cm, pp, G, sim)
                y = np.vstack([line_delta[tl][p] for p in pp]).reshape(-1)
                keep = np.broadcast_to(
                    (line_cm[tl] > args.min_expr) & line_cov[tl],
                    (len(pp), G)).copy()
                for j, p in enumerate(pp):
                    if p in gidx:
                        keep[j, gidx[p]] = False
                keep = keep.reshape(-1) & np.isfinite(y)
                X = X[keep].astype(np.float64)
                y = y[keep].astype(np.float64)
                if not len(y):
                    continue
                if not np.isfinite(X).all():
                    raise ValueError(f"Non-finite training features: {tl}")
                if sx is None:
                    sx = np.zeros(X.shape[1])
                    sy = 0.0
                    sxx = np.zeros((X.shape[1], X.shape[1]))
                    sxy = np.zeros(X.shape[1])
                wt = 1.0 / nr
                sx += wt * X.sum(0)
                sy += wt * y.sum()
                sxx += wt * (X.T @ X)
                sxy += wt * (X.T @ y)
                total += wt * len(y)
        if total == 0:
            print(f"Skipping {line}: no training rows", flush=True)
            continue
        mu, ym = sx / total, sy / total
        cov = sxx / total - np.outer(mu, mu)
        cov = (cov + cov.T) / 2
        sd = np.sqrt(np.maximum(np.diag(cov), 0))
        sd[sd < 1e-8] = 1.0
        pc = dict(mu=mu, sd=sd, ym=ym,
                  XtX=cov / np.outer(sd, sd),
                  Xty=(sxy / total - mu * ym) / sd, F=len(mu))

        # EVAL on the held-out line's test perts (proxy: cosine of predicted vs measured delta,
        # own target removed), restricted to covered+expressed genes
        test_perts = [p for p in folds[line]["test_perts"] if p in line_delta[line]]
        # Choose the panel without looking at measured effect sizes.
        if len(test_perts) > 300:
            test_perts = sorted(np.random.default_rng(args.seed).choice(
                test_perts, 300, replace=False).tolist())
        if len(test_perts) < 20:
            print(f"Skipping {line}: fewer than 20 test perts", flush=True)
            continue
        X, donors, fnames = build_rows(
            line, train_lines, line_delta, line_cov, line_cm, test_perts, G, sim)
        meas = np.vstack([line_delta[line][p] for p in test_perts])
        elig = (line_cm[line] > args.min_expr) & line_cov[line]
        for p in test_perts:
            if p in gidx:
                elig[gidx[p]] = False
        if elig.sum() < 2:
            print(f"Skipping {line}: no usable control profile/genes", flush=True)
            continue
        if not np.isfinite(X).all() or not np.isfinite(meas[:, elig]).all():
            raise ValueError(f"Non-finite evaluation data: {line}")

        def report(method, pred, **extra):
            P, M = pred[:, elig].astype(float), meas[:, elig].astype(float)
            def unit(A):
                return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-12)
            simmat = unit(P) @ unit(M).T
            correct = np.diag(simmat)[:, None]
            better = (simmat > correct + 1e-9).sum(1)
            ties = (np.abs(simmat - correct) <= 1e-9).sum(1) - 1
            pds = float(np.mean(1 - (better + 0.5 * ties) / (len(P) - 1)))
            rec = dict(
                line=line, seed=args.seed, method=method, n_test=len(P),
                median_cos=median_cos(pred, meas, test_perts, elig, gidx),
                pds_proxy=pds,
                mse_ratio=float(np.sum((P-M)**2) / max(np.sum(M*M), 1e-12)),
                **extra)
            rows.append(rec)
            print(rec, flush=True)

        report("no_change", np.zeros_like(meas))
        acc = np.zeros_like(meas, dtype=float)
        centered = np.zeros_like(acc)
        count = np.zeros_like(acc)
        for j, d in enumerate(donors):
            D = X[:, j].reshape(len(test_perts), G)
            present = np.array([p in line_delta[d] for p in test_perts])
            mask = present[:, None] & line_cov[d][None, :]
            acc += np.where(mask, D, 0)
            count += mask
            if present.any():
                center = D[present].mean(0)
                centered += np.where(mask, D - center, 0)
            report("single_" + d, D, coverage=float(present.mean()))
        for scale in (0.5, 0.7, 1.0):
            report(f"equal_raw_{scale}", scale * acc / np.maximum(count, 1))
            report(f"equal_centered_{scale}", scale * centered / np.maximum(count, 1))

        for alpha in alphas:
            w, b = ridge_solve(pc, alpha)
            pred = (X @ w + b).reshape(len(test_perts), G)
            report(f"ridge_{alpha}", pred)
            report(f"ridge_centered_{alpha}", pred - pred.mean(0))
            print("weights", dict(zip(fnames, np.round(w, 4))), flush=True)

    if args.out:
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
