#!/usr/bin/env python
"""Step 3b: find Challenge genes that PIE's gene_text source lacks, and write a filtered gene
file for `pie prep` (PIE errors on any gene axis entry without gene text).

Run in the `pie` conda env (needs the pie package; uses PIE's own loader so the key set is
exactly what `pie infer` checks against).

Usage:
    python scripts/06_pie_gene_filter.py --gene-names data_dir/gene_names.csv \
        --ckpt-config $PIE_RUNS_ROOT/replogle_xdataset/config.yaml \
        --out-genes data_dir/gene_names_pie.csv --out-missing data_dir/genes_missing_in_pie.txt

The filtered file has no header. The missing genes must later be filled with zero predicted
effect when assembling cells on the full 18,533-gene axis.
"""
import argparse

import yaml


def read_genes(path):
    genes = [line.strip().split(",")[0] for line in open(path) if line.strip()]
    if genes and genes[0].lower() in ("gene", "gene_name", "gene_names", "symbol"):
        genes = genes[1:]
    return genes


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gene-names", required=True)
    p.add_argument("--ckpt-config", required=True, help="config.yaml saved with the PIE checkpoint")
    p.add_argument("--out-genes", required=True)
    p.add_argument("--out-missing", required=True)
    args = p.parse_args()

    import pie.data.datamodule as dm  # read_source / resolve_asset are the loader pie infer uses

    cfg = yaml.safe_load(open(args.ckpt_config))
    gene_text_dir = cfg["data"]["gene_text_dir"]
    src = dm.read_source(dm.resolve_asset(gene_text_dir, kind="source"))

    genes = read_genes(args.gene_names)
    missing = [g for g in genes if g not in src.key_to_row]
    keep = [g for g in genes if g in src.key_to_row]
    print(f"{len(genes)} challenge genes; {len(missing)} missing from gene_text; {len(keep)} kept")
    print("first missing:", missing[:20])

    with open(args.out_genes, "w") as f:
        f.write("\n".join(keep) + "\n")
    with open(args.out_missing, "w") as f:
        f.write("\n".join(missing) + "\n")


if __name__ == "__main__":
    main()
