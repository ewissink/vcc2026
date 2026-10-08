#!/usr/bin/env python
"""Step 1a: download the scPerturb-prepared Replogle 2022 h5ad files from Zenodo
(record 13350497, v1.4), with MD5 verification. Standard library only.

The other datasets need manual download; see config/datasets.yaml.

Usage:
    python scripts/01_fetch_public.py --out /path/to/raw --which k562_essential rpe1
"""
import argparse
import hashlib
import urllib.request
from pathlib import Path

BASE = "https://zenodo.org/api/records/13350497/files/{}/content"
FILES = {
    "k562_essential": ("ReplogleWeissman2022_K562_essential.h5ad", "d8cba17576d1a8afc0f7d71b79cad0f7"),
    "k562_gwps": ("ReplogleWeissman2022_K562_gwps.h5ad", "13db594f8f1d2ccb88fec44a13e414dc"),
    "rpe1": ("ReplogleWeissman2022_rpe1.h5ad", "cc7f1ec50aeb3a3e1b4a6cfa713d80fa"),
}


def md5sum(path, chunk=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--which", nargs="+", default=["k562_essential", "rpe1"], choices=list(FILES))
    args = p.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for key in args.which:
        fname, md5 = FILES[key]
        dest = out / fname
        if dest.exists() and md5sum(dest) == md5:
            print(f"{fname}: already present, md5 ok")
            continue
        print(f"downloading {fname} ...")
        urllib.request.urlretrieve(BASE.format(fname), dest)
        got = md5sum(dest)
        if got != md5:
            raise SystemExit(f"{fname}: md5 mismatch (got {got}, expected {md5}); delete and retry")
        print(f"{fname}: md5 ok")


if __name__ == "__main__":
    main()
