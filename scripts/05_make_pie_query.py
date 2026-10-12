#!/usr/bin/env python
"""Step 3a: build the query file for `pie infer` from the Challenge perturbation list.

Format (pie AGENTS.md): {"<dataset>.<context>": [perturbation, ...]}. The dataset name and
context labels must match the controls-only dir written by `pie prep` (meta.json).

Usage:
    python scripts/05_make_pie_query.py --pert-counts data_dir/pert_counts.csv \
        --meta $PIE_DATA_ROOT/vcc2026_val/meta.json --out query.json [--n-perts 5]
"""
import argparse
import json

import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pert-counts", required=True)
    p.add_argument("--meta", required=True, help="meta.json from `pie prep controls_only=true`")
    p.add_argument("--out", required=True)
    p.add_argument("--n-perts", type=int, default=None, help="use only the first N (timing test)")
    args = p.parse_args()

    meta = json.load(open(args.meta))
    perts = pd.read_csv(args.pert_counts)["target_gene"].astype(str).tolist()
    if args.n_perts:
        perts = perts[: args.n_perts]
    query = {f"{meta['dataset']}.{ctx}": perts for ctx in meta["context_to_id"]}
    json.dump(query, open(args.out, "w"), indent=1)
    print(f"{len(query)} contexts x {len(perts)} perts -> {args.out}")
    print("keys:", list(query))


if __name__ == "__main__":
    main()
