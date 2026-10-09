#!/usr/bin/env python
"""Step 5a: build a local scoring harness from a public Perturb-seq screen with real cells.

From one Replogle-style h5ad (raw counts) it writes, into --out:
  gene_names.csv       genes kept (Challenge order, only genes present in the source)
  pert_counts.csv      the panel: --n-perts perturbations with the largest measured effects
  context_<NAME>.h5ad  half of the control cells; use as the 'context' for PIE and 09_make_cells.py
  real.h5ad            panel perturbed cells (<= --max-cells each) + the OTHER half of the controls
  ctrl_heldout.h5ad    the held-out controls alone (append to predictions before scoring)
Effect size comes from the 02_harmonize deltas (norm over covered genes, own target excluded), so
the panel resembles the Challenge's enriched-for-effect panels. The harness is only as
representative as its source: Replogle is 3'-capture, not 10x Flex.

Usage:
    python scripts/11_build_harness_data.py --h5ad raw/ReplogleWeissman2022_K562_essential.h5ad \
        --harmonized-npz results/harmonized/replogle_k562.npz --gene-names data_dir/gene_names.csv \
        --name K562 --out harness_K562 [--n-perts 300] [--max-cells 400] [--min-cells 50]
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_genes(path):
    genes = [line.strip().split(",")[0] for line in open(path) if line.strip()]
    if genes and genes[0].lower() in ("gene", "gene_name", "gene_names", "symbol"):
        genes = genes[1:]
    return genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--harmonized-npz", required=True)
    ap.add_argument("--gene-names", required=True)
    ap.add_argument("--name", required=True, help="context label, e.g. K562")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pert-col", default="gene")
    ap.add_argument("--control-label", default="non-targeting")
    ap.add_argument("--n-perts", type=int, default=300)
    ap.add_argument("--max-cells", type=int, default=400)
    ap.add_argument("--min-cells", type=int, default=50, help="panel perts need at least this many cells")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    z = np.load(args.harmonized_npz, allow_pickle=True)
    perts, delta, covered = z["perts"].astype(str), z["delta"], z["covered"]
    challenge = read_genes(args.gene_names)
    gidx = {g: i for i, g in enumerate(challenge)}

    a = ad.read_h5ad(args.h5ad, backed="r")
    src_genes = a.var_names.astype(str)
    keep_src = [j for j, g in enumerate(src_genes) if g in gidx]
    # Challenge order, restricted to genes in the source
    order = sorted(keep_src, key=lambda j: gidx[src_genes[j]])
    genes_kept = [src_genes[j] for j in order]
    cols = np.array(order)
    kept_idx = np.array([gidx[g] for g in genes_kept])
    print(f"{len(genes_kept)} genes kept of {len(challenge)} challenge genes", flush=True)

    obs_pert = a.obs[args.pert_col].astype(str).values
    counts = pd.Series(obs_pert).value_counts()
    eff = np.linalg.norm(delta[:, kept_idx], axis=1)           # includes own target; remove it:
    for i, p in enumerate(perts):
        if p in gidx and gidx[p] in set(kept_idx.tolist()):
            eff[i] = np.linalg.norm(np.delete(delta[i, kept_idx], np.where(kept_idx == gidx[p])[0]))
    ok = np.array([(counts.get(p, 0) >= args.min_cells) and (p in gidx) for p in perts])
    cand = np.where(ok)[0]
    panel = perts[cand[np.argsort(-eff[cand])][: args.n_perts]]
    print(f"panel: {len(panel)} perts, effect norm {eff[np.isin(perts, panel)].min():.2f}.."
          f"{eff[np.isin(perts, panel)].max():.2f} (all perts median {np.median(eff):.2f})", flush=True)

    ctrl_pos = np.where(obs_pert == args.control_label)[0]
    rng.shuffle(ctrl_pos)
    half = len(ctrl_pos) // 2
    ctx_pos, held_pos = np.sort(ctrl_pos[:half]), np.sort(ctrl_pos[half:])

    def take(pos):
        sub = a[pos].to_memory()
        X = sp.csr_matrix(sub.X)[:, cols].astype(np.int32)
        return X, sub.obs_names.astype(str)

    def make(X, labels, ctx, prefix):
        if isinstance(labels, str):
            labels = np.full(X.shape[0], labels)
        obs = pd.DataFrame({"target_gene": labels, "context": np.full(X.shape[0], ctx)})
        obs.index = [f"{prefix}{i}" for i in range(len(obs))]
        r = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=genes_kept))
        return r

    Xc, _ = take(ctx_pos)
    make(Xc, args.control_label, args.name, "ctx").write_h5ad(out / f"context_{args.name}.h5ad")
    Xh, _ = take(held_pos)
    held = make(Xh, args.control_label, args.name, "held")
    held.write_h5ad(out / "ctrl_heldout.h5ad")

    parts, labs = [Xh], [np.array([args.control_label] * Xh.shape[0])]
    for p in panel:
        pos = np.where(obs_pert == p)[0]
        if len(pos) > args.max_cells:
            pos = np.sort(rng.choice(pos, args.max_cells, replace=False))
        Xp, _ = take(pos)
        parts.append(Xp)
        labs.append(np.array([p] * Xp.shape[0]))
    real = make(sp.vstack(parts, format="csr"), np.concatenate(labs), args.name, "real")
    real.write_h5ad(out / "real.h5ad")

    pd.Series(genes_kept).to_csv(out / "gene_names.csv", header=False, index=False)
    pd.DataFrame({"target_gene": panel}).to_csv(out / "pert_counts.csv", index=False)
    print(f"context controls {Xc.shape[0]}, held-out controls {Xh.shape[0]}, "
          f"real {real.shape} -> {out}", flush=True)


if __name__ == "__main__":
    main()
