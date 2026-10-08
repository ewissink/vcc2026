#!/usr/bin/env python
"""Step 1a: fetch the public datasets that can be fetched programmatically.

Only the Replogle K562/RPE1 sets go through pertpy (scPerturb-prepared AnnData). The rest
need manual download; see config/datasets.yaml.

Usage:
    python scripts/01_fetch_public.py --out /path/to/raw --which k562_essential k562_gwps
"""
import argparse
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--which", nargs="+", default=["k562_essential"],
                   choices=["k562_essential", "k562_gwps", "rpe1"])
    args = p.parse_args()

    import pertpy as pt  # pip install pertpy

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    loaders = {
        "k562_essential": pt.data.replogle_2022_k562_essential,
        "k562_gwps": pt.data.replogle_2022_k562_gwps,
        "rpe1": pt.data.replogle_2022_rpe1,
    }
    for name in args.which:
        print(f"fetching {name} ...")
        adata = loaders[name]()
        adata.write_h5ad(out / f"replogle_{name}.h5ad")
        print(adata)
        print("obs columns:", list(adata.obs.columns))
        print("var columns:", list(adata.var.columns))


if __name__ == "__main__":
    main()
