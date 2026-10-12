#!/usr/bin/env python
"""Step 5c: score predicted-cell files against a real reference with cell-eval2's vcc2026 preset.

Prints one row per prediction file with the six scored metrics averaged over perturbations,
and writes per-perturbation values to --out-prefix_perpert.csv. Values are RAW (not rescaled to
the leaderboard's 0 = mean-response baseline / 1 = replicate), so compare settings with each
other and with a no-change reference run, not with leaderboard numbers. Direction:
higher is better for pds_cosine, direction_fidelity, direction_reach, sig_jaccard;
lower is better for expr_mse_unbiased_capped_norm and lfc_nmae.

Run in an env with `pip install cell-eval2` (Python >= 3.11; `pip install pdex` speeds up CPU DE).
Both files must contain 'non-targeting' control cells (14_append_controls.py).

Usage:
    python scripts/15_score_local.py --real harness_K562gw/real.h5ad \
        --pred nochange=pred_nochange_withctrl.h5ad s1.0=pred_h_s1.0_withctrl.h5ad \
               s0.5=pred_h_s0.5_withctrl.h5ad --out-prefix harness_scores
"""
import argparse
import dataclasses

import pandas as pd

SCORED = {
    "pds_cosine": "pds (higher)",
    "expr_mse_unbiased_capped_norm": "mse (lower)",
    "de_wilcoxon_direction_fidelity_yield_raw": "fid (higher)",
    "de_wilcoxon_direction_reach_raw": "reach (higher)",
    "de_wilcoxon_sig_jaccard": "jac (higher)",
    "de_wilcoxon_lfc_nmae": "nmae (lower)",
}


def to_pandas(df):
    return df.to_pandas() if hasattr(df, "to_pandas") else pd.DataFrame(df)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--pred", nargs="+", required=True, help="name=path pairs")
    ap.add_argument("--out-prefix", default="harness_scores")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--pert-col", default="target_gene")
    ap.add_argument("--control", default="non-targeting")
    ap.add_argument("--cache-real", default="cache_real", help="cell-eval2 cache dir for the reference")
    args = ap.parse_args()

    from cell_eval2 import EvalConfig, compute_metrics

    cfg = dataclasses.replace(EvalConfig.from_preset("vcc2026"), pert_col=args.pert_col,
                              control=args.control, device=args.device, cache_real=args.cache_real)
    rows, per_pert = [], []
    for item in args.pred:
        name, path = item.split("=", 1)
        print(f"scoring {name} ({path}) ...", flush=True)
        df = to_pandas(compute_metrics(path, args.real, config=cfg))
        df = df[df["metric"].isin(SCORED)].copy()
        df["setting"] = name
        per_pert.append(df)
        m = df.groupby("metric")["value"].mean()
        rows.append({"setting": name, **{SCORED[k]: m.get(k, float("nan")) for k in SCORED}})
    table = pd.DataFrame(rows).set_index("setting")
    pd.concat(per_pert).to_csv(f"{args.out_prefix}_perpert.csv", index=False)
    table.to_csv(f"{args.out_prefix}_summary.csv")
    pd.set_option("display.width", 200)
    print("\nmean over perturbations (raw values):")
    print(table.round(4).to_string())


if __name__ == "__main__":
    main()
