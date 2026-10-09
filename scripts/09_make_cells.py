#!/usr/bin/env python
"""Step 4: turn PIE predictions (p_de, lfc_pred per gene) into raw-count single cells.

For each (context, perturbation) row:
  1. sample --n-cells real control cells from that context (with replacement);
  2. effect e_g = lfc_scale * lfc_pred_g for genes with p_de_g > p-threshold, else 0;
  3. down (e<0): binomial thinning of each cell's counts with p = 2^e;
     up (e>0): add Poisson(mu_g * (2^e - 1) * cell depth factor), where mu_g is the context's
     mean control count; (multiplying counts cannot create expression from zero);
  4. genes PIE cannot predict (not on its axis) are left unchanged.
Output: one sparse int32 h5ad with obs [context, target_gene], 18,533 genes in gene_names order,
checked against the Challenge limits before writing. p-threshold and lfc-scale are placeholders
that need calibration.

Usage:
    python scripts/09_make_cells.py --pred pred_val_full.parquet --pie-meta $PIE_DATA_ROOT/vcc2026_val/meta.json \
        --data-dir data_dir --pert-counts data_dir/pert_counts.csv --out prediction.h5ad \
        [--p-threshold 0.5] [--lfc-scale 1.0] [--workers 32] [--seed 0] [--n-cells 400]
"""
import argparse
import json
import multiprocessing as mp
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.sparse as sp

G = {}  # shared with forked workers


def read_genes(path):
    genes = [line.strip().split(",")[0] for line in open(path) if line.strip()]
    if genes and genes[0].lower() in ("gene", "gene_name", "gene_names", "symbol"):
        genes = genes[1:]
    return genes


def make_pop(job):
    ctx, pert, p_de, lfc, seed = job
    X, mu, depth = G["X"][ctx], G["mu"][ctx], G["depth"][ctx]
    rng = np.random.default_rng(seed)
    n_genes = X.shape[1]
    e = np.zeros(n_genes)
    sel = p_de > G["p_thr"]
    e[G["pie_to_full"][sel]] = G["lfc_scale"] * lfc[sel]
    f = 2.0 ** e
    idx = rng.integers(0, X.shape[0], size=G["n_cells"])
    sub = X[idx].tocsr()
    sub.data = rng.binomial(sub.data.astype(np.int64), np.minimum(f, 1.0)[sub.indices]).astype(np.int32)
    up = np.where(f > 1.0)[0]
    if len(up):
        lam = mu[up] * (f[up] - 1.0) * (depth[idx] / depth.mean())[:, None]
        add = rng.poisson(lam).astype(np.int32)
        r, c = np.nonzero(add)
        sub = sub + sp.csr_matrix((add[r, c], (r, up[c])), shape=sub.shape, dtype=np.int32)
    sub.eliminate_zeros()
    return ctx, pert, sub.astype(np.int32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--pie-meta", required=True, help="meta.json of the controls-only dir (gene axis)")
    ap.add_argument("--data-dir", required=True, help="dir with context_*.h5ad and gene_names.csv")
    ap.add_argument("--pert-counts", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--p-threshold", type=float, default=0.5)
    ap.add_argument("--lfc-scale", type=float, default=1.0)
    ap.add_argument("--n-cells", type=int, default=400)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    genes = read_genes(Path(args.data_dir) / "gene_names.csv")
    gidx = {g: i for i, g in enumerate(genes)}
    pie_genes = json.load(open(args.pie_meta))["genes"]
    G["pie_to_full"] = np.array([gidx[g] for g in pie_genes])
    pred = pq.read_table(args.pred).to_pandas()
    pert_order = pd.read_csv(args.pert_counts)["target_gene"].astype(str).tolist()
    G["p_thr"], G["lfc_scale"], G["n_cells"] = args.p_threshold, args.lfc_scale, args.n_cells
    print(f"{len(pred)} predictions; PIE axis {len(pie_genes)} / {len(genes)} genes", flush=True)

    G["X"], G["mu"], G["depth"] = {}, {}, {}
    for ctx in sorted(pred["context"].unique()):
        a = ad.read_h5ad(Path(args.data_dir) / f"context_{ctx}.h5ad")
        assert list(a.var_names) == genes, f"gene order mismatch in context {ctx}"
        X = sp.csr_matrix(a.X).astype(np.int32)
        G["X"][ctx] = X
        G["mu"][ctx] = np.asarray(X.mean(axis=0)).ravel()
        G["depth"][ctx] = np.asarray(X.sum(axis=1)).ravel().astype(float)

    lookup = {(r.context, r.perturbation): (np.asarray(r.p_de), np.asarray(r.lfc_pred))
              for r in pred.itertuples()}
    jobs, k = [], 0
    for ctx in sorted(G["X"]):
        for pert in pert_order:
            p_de, lfc = lookup[(ctx, pert)]
            assert len(p_de) == len(pie_genes), "prediction length != PIE gene axis"
            jobs.append((ctx, pert, p_de, lfc, [args.seed, k]))
            k += 1

    with mp.get_context("fork").Pool(args.workers) as pool:
        results = pool.map(make_pop, jobs, chunksize=4)

    X = sp.vstack([r[2] for r in results], format="csr")
    obs = pd.DataFrame({"context": np.repeat([r[0] for r in results], args.n_cells),
                        "target_gene": np.repeat([r[1] for r in results], args.n_cells)})
    row_sums = np.asarray(X.sum(axis=1)).ravel()
    assert X.shape[1] == len(genes) and X.dtype == np.int32 and X.data.min() >= 0
    assert row_sums.max() <= 1_000_000, "a cell exceeds 1,000,000 counts"
    assert X.nnz <= 4_750_000_000, "too many stored entries"
    assert (obs["target_gene"] != "non-targeting").all()
    print(f"cells {X.shape}, stored entries {X.nnz} ({X.nnz / X.shape[0]:.0f}/cell), "
          f"median UMI {np.median(row_sums):.0f}", flush=True)
    ad.AnnData(X=X, obs=obs.set_index(np.arange(len(obs)).astype(str)),
               var=pd.DataFrame(index=genes)).write_h5ad(args.out)
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
