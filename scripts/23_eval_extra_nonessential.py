"""Evaluate extra donors on held-out genome-wide panels.
Essential-screen membership is an exclusion proxy.
Metrics and candidate weights are exploratory, not official scores.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

ESSENTIAL = [
    "replogle_k562", "replogle_rpe1",
    "nadig_hepg2", "nadig_jurkat",
]
BASE = [
    "replogle_k562_gwps", "xatlas_hct116",
    "xatlas_hek293t", "cd4_rest",
]
EXTRA = ["pan_h1", "feng_sf", "feng_snf"]
TARGETS = BASE[:3]


def load(path, wanted):
    with np.load(path, allow_pickle=True) as z:
        perts = z["perts"].astype(str)
        if len(set(perts)) != len(perts):
            raise ValueError(f"Duplicate perturbations: {path}")
        hit = [i for i, p in enumerate(perts) if p in wanted]
        full = z["delta"]
        delta = full[hit].astype(np.float32)
        del full
        cov = z["covered"].astype(bool)
        if not np.isfinite(delta[:, cov]).all():
            raise ValueError(f"Non-finite covered effects: {path}")
        return dict(
            delta=delta, cov=cov, genes=z["genes"].astype(str),
            index={perts[j]: i for i, j in enumerate(hit)},
        )


def panel_rows(d, panel, center):
    A = np.zeros((len(panel), len(d["genes"])), dtype=np.float64)
    hit = np.array([p in d["index"] for p in panel])
    for i, p in enumerate(panel):
        if hit[i]:
            A[i] = d["delta"][d["index"][p]]
    if center and hit.any():
        A[hit] -= A[hit].mean(axis=0)
    mask = hit[:, None] & d["cov"][None, :]
    return np.where(mask, A, 0), mask, hit


def metrics(P, M):
    def unit(A):
        norm = np.linalg.norm(A, axis=1, keepdims=True)
        return A / np.maximum(norm, 1e-12)

    S = unit(P) @ unit(M).T
    true = np.diag(S)[:, None]
    better = (S > true + 1e-9).sum(axis=1)
    ties = (np.abs(S - true) <= 1e-9).sum(axis=1) - 1
    pds = 1 - (better + .5 * ties) / (len(P) - 1)
    error = ((P - M) ** 2).sum(axis=1)
    truth = (M ** 2).sum(axis=1)
    return pds, error, truth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    folder = args.root / "results/harmonized"

    with (args.root / "data_dir/pert_counts.csv").open() as f:
        official = set(r["target_gene"].strip() for r in csv.DictReader(f))

    excluded = set()
    for name in ESSENTIAL:
        with np.load(folder / f"{name}.npz", allow_pickle=True) as z:
            excluded.update(z["perts"].astype(str))

    print(f"Essential-screen exclusion set: {len(excluded)} genes", flush=True)
    print(
        f"Competition targets excluded: "
        f"{len(official & excluded)}/{len(official)}",
        flush=True,
    )
    for name in EXTRA:
        with np.load(folder / f"{name}.npz", allow_pickle=True) as z:
            p = set(z["perts"].astype(str))
        print(
            f"{name}: competition overlap before={len(p & official)}, "
            f"after={len(p & (official - excluded))}; "
            f"screen targets remaining={len(p - excluded)}/{len(p)}",
            flush=True,
        )

    folds = json.loads((args.root / "results/splits.json").read_text())
    panels = {}
    for name in TARGETS:
        with np.load(folder / f"{name}.npz", allow_pickle=True) as z:
            line = str(z["cell_line"])
            genes = z["genes"].astype(str)
            gidx = {g: i for i, g in enumerate(genes)}
            cm = z["ctrl_mean"]
            cov = z["covered"].astype(bool)
            pool = sorted(
                (set(z["perts"].astype(str))
                 & set(folds[line]["test_perts"])) - excluded
            )
            pool = [
                p for p in pool
                if p in gidx and cov[gidx[p]]
                and np.isfinite(cm[gidx[p]]) and cm[gidx[p]] > .5
            ]

        if len(pool) < 20:
            raise ValueError(f"{name}: only {len(pool)} eligible test targets")
        panel = sorted(np.random.default_rng(args.seed).choice(
            pool, min(300, len(pool)), replace=False
        ).tolist())
        panels[name] = (
            panel, np.isfinite(cm) & (cm > .5) & cov
        )
        print(
            f"{name}: eligible test pool={len(pool)}, panel={len(panel)}",
            flush=True,
        )

    # Retain only selected rows; decompress one source matrix at a time.
    wanted = set().union(*(set(p) for p, _ in panels.values()))
    donors = {}
    for name in BASE + EXTRA:
        print(f"Loading selected rows: {name}", flush=True)
        donors[name] = load(folder / f"{name}.npz", wanted)

    genes = donors[BASE[0]]["genes"]
    for name, d in donors.items():
        if not np.array_equal(genes, d["genes"]):
            raise ValueError(f"Gene-axis mismatch: {name}")

    records = []
    for target, (panel, expressed) in panels.items():
        eligible = expressed & ~np.isin(genes, panel)
        if eligible.sum() < 100:
            raise ValueError(f"{target}: too few eligible genes")

        M = panel_rows(donors[target], panel, False)[0][:, eligible]
        total = np.zeros_like(M)
        count = np.zeros_like(M)

        for name in BASE:
            if name == target:
                continue
            A, mask, _ = panel_rows(donors[name], panel, True)
            total += A[:, eligible]
            count += mask[:, eligible]

        baseline = .7 * total / np.maximum(count, 1)
        bp, be, _ = metrics(baseline, M)

        for extra in EXTRA:
            A, mask, hit = panel_rows(donors[extra], panel, True)
            A, mask = A[:, eligible], mask[:, eligible]
            denominator = np.maximum(count + mask, 1)

            predictions = {
                "baseline": baseline,
                "shrink_only": .7 * total / denominator,
            }
            for beta in (.25, 1., 4.):
                predictions[f"add_beta_{beta}"] = (
                    .7 * (total + beta * A) / denominator
                )

            for method, P in predictions.items():
                pds, error, truth = metrics(P, M)
                subsets = [
                    ("all", np.ones(len(panel), bool)),
                    ("extra_covered", hit),
                ]
                for subset, keep in subsets:
                    n = int(keep.sum())
                    if not n:
                        continue
                    den = float(truth[keep].sum())
                    rec = dict(
                        target=target,
                        seed=args.seed,
                        extra=extra,
                        method=method,
                        subset=subset,
                        n=n,
                        n_panel=len(panel),
                        pds_proxy=float(pds[keep].mean()),
                        pds_change_vs_baseline=float(
                            (pds[keep] - bp[keep]).mean()
                        ),
                        mse_ratio=(
                            float(error[keep].sum() / den)
                            if den > 0 else float("nan")
                        ),
                        mse_change_vs_baseline=(
                            float((error[keep] - be[keep]).sum() / den)
                            if den > 0 else float("nan")
                        ),
                    )
                    records.append(rec)
                    print(rec, flush=True)

    out = (
        args.root
        / f"results/transfer_fixed/nonessential_seed{args.seed}.csv"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(f"Wrote {out}", flush=True)


if __name__ == "__main__":
    main()
