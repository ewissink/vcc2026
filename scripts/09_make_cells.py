#!/usr/bin/env python
"""Step 4: turn PIE predictions (p_de, lfc_pred per gene) into raw-count single cells.

For each (context, perturbation) row:
  1. sample --n-cells real control cells from that context (with replacement);
  2. effect e_g = lfc_scale * clip(lfc_pred_g, +-lfc-clip) for genes with p_de_g > p-threshold AND
     mean control CPM in that context > min-cpm (the metric only tests such genes, and PIE's
     fold changes on near-unexpressed genes are huge artifacts), else 0;
  3. optionally blend toward the context's mean effect (--shrink) and cap the genes changed
     per perturbation (--max-genes);
  4. down (e<0): binomial thinning of each cell's counts with p = 2^e;
     up (e>0): add Poisson(mu_g * (2^e - 1) * cell depth factor), where mu_g is the context's
     mean control count; (multiplying counts cannot create expression from zero);
  4b. optionally add a donor screen's measured effects (--donor-npz from 02_harmonize.py):
     e += donor-scale * (delta / ln2), restricted to genes the donor covers and that are expressed
     in the context; delta is centered on the donor's mean over the panel (--no-donor-center to
     keep it) so the donor's own shared response is not imported; perturbations absent from the
     donor are left at the PIE base; the total is clipped to +-lfc-clip.
  5. genes PIE cannot predict (not on its axis) are left unchanged.
Output: one sparse int32 h5ad with obs [context, target_gene], 18,533 genes in gene_names order,
checked against the Challenge limits before writing. p-threshold and lfc-scale are placeholders
that need calibration.

Usage:
    python scripts/09_make_cells.py --pred pred_val_full.parquet --pie-meta $PIE_DATA_ROOT/vcc2026_val/meta.json \
        --data-dir data_dir --pert-counts data_dir/pert_counts.csv --out prediction.h5ad \
        [--p-threshold 0.5] [--lfc-scale 1.0] [--min-cpm 5] [--lfc-clip 3.0]
        [--shrink 1.0] [--max-genes N]
        [--donor-npz results/harmonized/replogle_k562_gwps.npz --donor-scale 0.1] [--workers 32] [--seed 0] [--n-cells 400]
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
    ctx, pert, e, seed = job
    X, mu, depth = G["X"][ctx], G["mu"][ctx], G["depth"][ctx]
    rng = np.random.default_rng(seed)
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
    ap.add_argument("--shrink", type=float, default=1.0,
                    help="blend toward the context's mean effect: 1 = PIE per-pert, 0 = mean only")
    ap.add_argument("--max-genes", type=int, default=None,
                    help="after blending, keep only the N strongest genes per perturbation")
    ap.add_argument("--donor-npz", default=None, help="harmonized donor screen to blend in")
    ap.add_argument("--donor-scale", type=float, default=0.1)
    ap.add_argument("--no-donor-center", action="store_true")
    ap.add_argument("--donor-contexts", nargs="+", default=None,
                    help="apply the donor only to these contexts (default: all)")
    ap.add_argument("--min-cpm", type=float, default=5.0)
    ap.add_argument("--lfc-clip", type=float, default=3.0, help="clip |log2 fc| (log2 units)")
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
    G["n_cells"] = args.n_cells
    print(f"{len(pred)} predictions; PIE axis {len(pie_genes)} / {len(genes)} genes", flush=True)

    G["X"], G["mu"], G["depth"], G["expr_ok"], G["expr_full"] = {}, {}, {}, {}, {}
    for ctx in sorted(pred["context"].unique()):
        a = ad.read_h5ad(Path(args.data_dir) / f"context_{ctx}.h5ad")
        assert list(a.var_names) == genes, f"gene order mismatch in context {ctx}"
        X = sp.csr_matrix(a.X).astype(np.int32)
        G["X"][ctx] = X
        G["mu"][ctx] = np.asarray(X.mean(axis=0)).ravel()
        G["depth"][ctx] = np.asarray(X.sum(axis=1)).ravel().astype(float)
        cpm = np.asarray((sp.diags(1e6 / np.maximum(G["depth"][ctx], 1)) @ X).mean(axis=0)).ravel()
        G["expr_full"][ctx] = cpm > args.min_cpm
        G["expr_ok"][ctx] = G["expr_full"][ctx][G["pie_to_full"]]
        print(f"context {ctx}: {int(G['expr_ok'][ctx].sum())} PIE genes pass min-cpm {args.min_cpm}",
              flush=True)

    lookup = {(r.context, r.perturbation): (np.asarray(r.p_de), np.asarray(r.lfc_pred))
              for r in pred.itertuples()}
    nsel = np.array([int(((np.asarray(r.p_de) > args.p_threshold) & G["expr_ok"][r.context]).sum())
                     for r in pred.itertuples()])
    print(f"genes changed per perturbation after filters: median {int(np.median(nsel))}, "
          f"90th pct {int(np.percentile(nsel, 90))}, max {int(nsel.max())}", flush=True)
    def effect(ctx, p_de, lfc):
        e = np.zeros(len(genes))
        sel = (p_de > args.p_threshold) & G["expr_ok"][ctx]
        e[G["pie_to_full"][sel]] = args.lfc_scale * np.clip(lfc[sel], -args.lfc_clip, args.lfc_clip)
        return e

    donor = None
    if args.donor_npz:
        z = np.load(args.donor_npz, allow_pickle=True)
        assert list(z["genes"].astype(str)) == genes, "donor gene axis != Challenge gene order"
        d_perts = z["perts"].astype(str)
        d_idx = {g: i for i, g in enumerate(d_perts)}
        d_delta = z["delta"].astype(np.float64) / np.log(2.0)          # log1p-CPM units -> log2
        d_cov = z["covered"].astype(bool)
        panel_rows = [d_idx[g] for g in pert_order if g in d_idx]
        d_mean = d_delta[panel_rows].mean(axis=0) if panel_rows else np.zeros(len(genes))
        print(f"donor {Path(args.donor_npz).stem}: {len(panel_rows)}/{len(pert_order)} panel perts "
              f"covered, scale {args.donor_scale}, centered={not args.no_donor_center}", flush=True)
        donor = True

    jobs, k = [], 0
    for ctx in sorted(G["X"]):
        E = np.vstack([effect(ctx, *lookup[(ctx, pert)]) for pert in pert_order])
        assert all(len(lookup[(ctx, p)][0]) == len(pie_genes) for p in pert_order)
        mean_e = E.mean(axis=0)                       # the context's shared response
        E = mean_e + args.shrink * (E - mean_e)       # shrink=1 keeps PIE; 0 = mean only
        if args.max_genes:
            for i in range(E.shape[0]):               # keep the strongest N per perturbation
                drop = np.argsort(-np.abs(E[i]))[args.max_genes:]
                E[i, drop] = 0.0
        if donor and (args.donor_contexts is None or ctx in args.donor_contexts):
            ok = d_cov & G["expr_full"][ctx]
            for i, pert in enumerate(pert_order):
                if pert in d_idx:
                    dev = d_delta[d_idx[pert]] - (0.0 if args.no_donor_center else d_mean)
                    E[i, ok] += args.donor_scale * dev[ok]
            E = np.clip(E, -args.lfc_clip, args.lfc_clip)
        print(f"context {ctx}: genes changed per pert after shrink/cap: "
              f"median {int(np.median((E != 0).sum(1)))}", flush=True)
        for i, pert in enumerate(pert_order):
            jobs.append((ctx, pert, E[i], [args.seed, k]))
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
