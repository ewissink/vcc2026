#!/usr/bin/env python
"""Step 1e: shrink each Challenge perturbation's donor delta by its estimated reliability.

Most K562 GWPS deltas for the Challenge perts are sampling noise (split-half check: ~100/272
perts beat a control-cell null). Writing delta = signal + noise, the MSE-optimal shrinkage of the
delta toward 0 is f = S / (S + N), the signal fraction of its squared norm. Two estimates of f,
from the split-half output (the in-memory split-half script; columns set, pert, cells, eff,
split_cos, top200_sign; set = challenge | null):
  cos:   split-half cosine c ~ S / (S + 2N)  ->  f = 2c / (1 + c)   (Spearman-Brown), c < 0 -> 0
  norm:  noise norm^2 of a full delta ~ K / n_cells, K = median(eff^2 * cells) over the null
         control groups  ->  f = 1 - K / (n_cells * eff^2)
  mean:  average of the two (default; cos is noisy per pert, norm assumes perturbed cells are
         as noisy as controls).
Each listed perturbation's delta row is multiplied by its f (own target gene kept at its
measured value unless --no-keep-target); perturbations not in the CSV are left unchanged (09
only uses panel perts). The donor transfer scale is still set in 09 (--donor-scale / --donor).

Usage:
    python scripts/17_reliability_weight.py --npz results/harmonized/replogle_k562_gwps.npz \
        --split-csv split_half_challenge.csv --out results/weighted/replogle_k562_gwps_rw.npz \
        [--method mean|cos|norm] [--min-weight 0]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def weights(rows, K, method):
    c = rows["split_cos"].clip(lower=0.0)
    f_cos = (2 * c / (1 + c)).clip(0.0, 1.0)
    f_norm = (1 - K / (rows["cells"] * rows["eff"] ** 2)).clip(0.0, 1.0)
    return {"cos": f_cos, "norm": f_norm, "mean": (f_cos + f_norm) / 2}[method], f_cos, f_norm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True, help="harmonized donor from 02_harmonize.py")
    ap.add_argument("--split-csv", required=True, help="split_half_challenge.csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--method", choices=["mean", "cos", "norm"], default="mean")
    ap.add_argument("--min-weight", type=float, default=0.0, help="floor on f")
    ap.add_argument("--no-keep-target", action="store_true")
    args = ap.parse_args()

    # keep_default_na=False: pandas would otherwise read the label "null" (and gene names such as
    # "NA") as missing values
    df = pd.read_csv(args.split_csv, keep_default_na=False, na_values=["", "nan", "NaN"])
    ch = df[df["set"] == "challenge"].drop_duplicates("pert").set_index("pert")
    null = df[df["set"] == "null"]
    if ch.empty or null.empty:
        raise SystemExit("split CSV needs both 'challenge' and 'null' rows")
    K = float(np.median(null["eff"] ** 2 * null["cells"]))
    f, f_cos, f_norm = weights(ch, K, args.method)
    f = f.clip(lower=args.min_weight)
    fn, _, _ = weights(null, K, args.method)            # what pure-noise groups would get

    q = [0.1, 0.25, 0.5, 0.75, 0.9]
    print(f"noise constant K = {K:.4g} (median eff^2 * cells over {len(null)} null groups)")
    print(f"agreement of the two estimates: Pearson(cos, norm) = {np.corrcoef(f_cos, f_norm)[0, 1]:.3f}")
    print(f"weights ({args.method}) over {len(f)} Challenge perts, quantiles {q}: "
          f"{np.round(np.quantile(f, q), 3).tolist()}")
    print(f"  perts with f > 0.25: {int((f > 0.25).sum())}, f > 0.5: {int((f > 0.5).sum())}, "
          f"f == 0: {int((f == 0).sum())}")
    print(f"null groups (should be ~0): median {np.median(fn):.3f}, 90th pct {np.quantile(fn, 0.9):.3f}")

    z = np.load(args.npz, allow_pickle=True)
    data = {k: z[k] for k in z.files}
    perts = data["perts"].astype(str)
    genes = data["genes"].astype(str)
    cov = data["covered"].astype(bool)
    delta = data["delta"].astype(np.float32).copy()
    pidx = {p: i for i, p in enumerate(perts)}
    gidx = {g: i for i, g in enumerate(genes)}
    done = []
    for p, w in f.items():
        i = pidx.get(p)
        if i is None:
            continue
        t = gidx.get(p)
        keep = t is not None and cov[t] and not args.no_keep_target
        tv = delta[i, t] if keep else None
        delta[i] *= np.float32(w)
        if keep:
            delta[i, t] = tv
        done.append((p, float(w)))
    print(f"weighted {len(done)}/{len(f)} Challenge perts present in the donor "
          f"({len(perts) - len(done)} other donor perts unchanged)")

    data["delta"] = delta
    data["reliability_perts"] = np.array([p for p, _ in done])
    data["reliability_weights"] = np.array([w for _, w in done], dtype=np.float32)
    data["reliability_method"] = args.method
    data["reliability_source"] = str(args.split_csv)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, **data)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
