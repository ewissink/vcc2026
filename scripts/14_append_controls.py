#!/usr/bin/env python
"""Step 5b: append held-out control cells to a predicted-cells h5ad so cell-eval2 can score it.

The Challenge platform adds real held-out controls to submissions itself; locally cell-eval2
needs 'non-targeting' cells in both the predicted and the reference file (see 10_celleval_smoke.py).

Usage:
    python scripts/14_append_controls.py --pred pred_h_s0.5.h5ad \
        --controls harness_K562gw/ctrl_heldout.h5ad --out pred_h_s0.5_withctrl.h5ad
"""
import argparse

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--controls", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pred = ad.read_h5ad(args.pred)
    ctrl = ad.read_h5ad(args.controls)
    assert list(pred.var_names) == list(ctrl.var_names), "gene axes differ between pred and controls"
    assert (ctrl.obs["target_gene"] == "non-targeting").all(), "controls file has non-control cells"
    assert not (pred.obs["target_gene"] == "non-targeting").any(), "pred already contains controls"

    X = sp.vstack([sp.csr_matrix(pred.X), sp.csr_matrix(ctrl.X)], format="csr").astype(np.int32)
    obs = pd.concat([pred.obs[["target_gene", "context"]], ctrl.obs[["target_gene", "context"]]])
    obs.index = [f"c{i}" for i in range(len(obs))]
    ad.AnnData(X=X, obs=obs, var=pred.var).write_h5ad(args.out)
    print(f"{pred.n_obs} predicted + {ctrl.n_obs} control cells -> {args.out}")


if __name__ == "__main__":
    main()
