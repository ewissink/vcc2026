#!/usr/bin/env python
"""Step 1c: leave-one-cell-line-out folds over harmonized datasets.

For each held-out cell line L (datasets with cell_line == L):
  - train_datasets: every dataset from other cell lines
  - eval perturbations of L are split into
      "seen_elsewhere": the gene is also perturbed in some training line (cross-context)
      "novel":          the gene is perturbed in no training line (zero-shot gene)
  - a random `--val-frac` of eval perts (stratified by the above) is the tuning set for
    ensemble weights; the rest is the test set for the fold.

Usage:
    python scripts/03_make_splits.py --harmonized harmonized/ --out splits.json [--seed 0]
"""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--harmonized", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--val-frac", type=float, default=0.3)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    ds = {}
    for f in sorted(Path(args.harmonized).glob("*.npz")):
        z = np.load(f, allow_pickle=False)
        ds[f.stem] = dict(cell_line=str(z["cell_line"]), perts=set(z["perts"].tolist()))
    lines = sorted({d["cell_line"] for d in ds.values()})
    print("datasets:", {k: (v["cell_line"], len(v["perts"])) for k, v in ds.items()})

    folds = {}
    for line in lines:
        held = [k for k, v in ds.items() if v["cell_line"] == line]
        train = [k for k in ds if k not in held]
        train_perts = set().union(*[ds[k]["perts"] for k in train]) if train else set()
        eval_perts = sorted(set().union(*[ds[k]["perts"] for k in held]))
        seen = [g for g in eval_perts if g in train_perts]
        novel = [g for g in eval_perts if g not in train_perts]
        val, test = [], []
        for grp in (seen, novel):
            grp = list(grp)
            rng.shuffle(grp)
            k = int(round(args.val_frac * len(grp)))
            val += grp[:k]
            test += grp[k:]
        folds[line] = dict(
            held_out_datasets=held, train_datasets=train,
            n_seen_elsewhere=len(seen), n_novel=len(novel),
            val_perts=sorted(val), test_perts=sorted(test),
        )
        print(f"{line}: train={train} eval={len(eval_perts)} seen={len(seen)} novel={len(novel)}")
    Path(args.out).write_text(json.dumps(folds, indent=1))


if __name__ == "__main__":
    main()
