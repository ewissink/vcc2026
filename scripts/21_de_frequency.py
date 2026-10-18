#!/usr/bin/env python
"""Step 1h: build a DE-frequency shared-response profile from donor screens (Liu/Wang 2026).

Their 3rd-place VCC-2025 recipe = DE-frequency ranking + cross-dataset transfer. DE frequency is a
perturbation-AGNOSTIC hierarchy: rank each outcome gene by how often it is DE across perturbations.
It is among the strongest predictors of which genes respond for a held-out perturbation.

We turn that into a shared-response vector to use in 09_make_cells (--shared-npz) in place of PIE's
panel-mean shared part (which has ~coin-flip directions). For each Challenge gene g, pooled over all
perturbations of all given donors:
  freq_g = fraction of perturbations in which g is DE
  mean_g = mean signed effect of g over perturbations where it is DE (its typical direction/size)
  shared_g = freq_g * mean_g                      (down-weights rarely-responsive genes)
DE is called from effect_se when present (|effect/se| >= --z), else from |effect| >= --delta-cut.
Units match the donors (log1p-CP delta); 09's --shared-npz divides by ln2 like any donor.

Use the high-coverage donors most like the Challenge contexts (K562/HCT116/HEK293T/CD4).

Usage:
    python scripts/21_de_frequency.py --donors results/harmonized/replogle_k562_gwps.npz \
        results/harmonized/xatlas_hct116.npz results/harmonized/xatlas_hek293t.npz \
        results/harmonized/cd4_rest.npz --gene-names data_dir/gene_names.csv \
        --out results/harmonized/de_freq_shared.npz
"""
import argparse
from pathlib import Path

import numpy as np


def read_genes(path):
    genes = [l.strip().split(",")[0] for l in open(path) if l.strip()]
    if genes and genes[0].lower() in ("gene", "gene_name", "gene_names", "symbol"):
        genes = genes[1:]
    return genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--gene-names", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--z", type=float, default=3.0, help="DE if |effect/se| >= z (se donors)")
    ap.add_argument("--delta-cut", type=float, default=0.5,
                    help="DE if |effect| >= this (donors without se)")
    ap.add_argument("--freq-power", type=float, default=1.0,
                    help="weight genes by freq**power (0 = plain mean effect, 1 = freq-weighted)")
    args = ap.parse_args()

    challenge = read_genes(args.gene_names)
    G = len(challenge)
    n_de = np.zeros(G)            # perts where gene is DE
    n_meas = np.zeros(G)          # perts where gene is covered
    sum_de = np.zeros(G)          # summed signed effect over DE perts

    for f in args.donors:
        z = np.load(f, allow_pickle=True)
        assert list(z["genes"].astype(str)) == challenge, f"{f}: gene axis != Challenge order"
        delta = z["delta"].astype(np.float64)                       # perts x G, log1p-CP delta
        cov = z["covered"].astype(bool)
        if "effect_se" in z.files:
            se = z["effect_se"].astype(np.float64)
            de = (np.abs(delta) >= args.z * np.where(se > 0, se, np.inf)) & cov[None, :]
        else:
            de = (np.abs(delta) >= args.delta_cut) & cov[None, :]
        n_meas += cov * delta.shape[0]
        n_de += de.sum(axis=0)
        sum_de += np.where(de, delta, 0.0).sum(axis=0)
        print(f"{Path(f).stem}: {delta.shape[0]} perts, median DE/pert "
              f"{int(np.median(de.sum(axis=1)))}", flush=True)

    freq = np.divide(n_de, n_meas, out=np.zeros(G), where=n_meas > 0)
    mean_de = np.divide(sum_de, n_de, out=np.zeros(G), where=n_de > 0)   # mean signed effect when DE
    shared = (freq ** args.freq_power) * mean_de
    covered = n_meas > 0

    q = [0.5, 0.9, 0.99]
    print(f"genes with freq>0: {int((freq>0).sum())}/{G}; freq quantiles "
          f"{np.round(np.quantile(freq[covered], q), 3).tolist()} at {q}", flush=True)
    print(f"shared profile nonzero {int((shared!=0).sum())}, "
          f"L2 norm {np.linalg.norm(shared):.2f}", flush=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, genes=np.array(challenge, dtype=str),
                        shared=shared.astype(np.float32), freq=freq.astype(np.float32),
                        covered=covered)
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
