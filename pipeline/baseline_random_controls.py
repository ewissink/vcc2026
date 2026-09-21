"""
Baseline submission -- random control resampling
==================================================

For each (cell_type, perturbation) pair required by the challenge,
draws 400 cells at random (with replacement if needed) from that
cell type's own real control-cell pool, and submits them unmodified
as the "predicted" perturbed population.

This makes no attempt to model any perturbation effect -- it is the
"predict no change" floor. It's included here because it's a
meaningful floor to know, not a target: per the blog analysis
discussed earlier, this kind of resampled-real-cells baseline
outperforms a naive identical-mean-profile submission (which scores
catastrophically low due to zero variance), but any actual modeling
(GEARS, STATE, Stack, or the ensemble) should be validated against
this floor to confirm it's adding real signal rather than just
matching or underperforming it.

Usage:
    python baseline_random_controls.py \
        --controls-dir /path/to/context_controls/ \
        --targets-csv /path/to/required_targets.csv \
        --out submission_baseline.h5ad \
        --n-cells 400 \
        --seed 0

Expected inputs:
    --controls-dir: one .h5ad per cell type, named <cell_type>.h5ad,
        each containing only that context's control (non-targeting)
        cells, raw counts in .X, genes in .var_names.
    --targets-csv: two columns, cell_type,perturbation -- the full
        list of (context, gene) pairs the submission must cover.
        This is the manifest the challenge provides; adjust the
        column names below if yours differ.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def load_control_pools(controls_dir: Path) -> dict[str, sc.AnnData]:
    """One AnnData of control cells per cell type, keyed by filename stem."""
    pools = {}
    for h5ad_path in sorted(controls_dir.glob("*.h5ad")):
        cell_type = h5ad_path.stem
        pools[cell_type] = sc.read_h5ad(h5ad_path)
    if not pools:
        raise FileNotFoundError(f"No .h5ad files found in {controls_dir}")
    return pools


def sample_baseline_population(
    control_adata: sc.AnnData,
    n_cells: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Randomly draws n_cells real control cells, with replacement if the
    control pool is smaller than n_cells. Returns raw counts, dtype
    matching the source (should already be integer counts).
    """
    n_available = control_adata.n_obs
    replace = n_available < n_cells
    idx = rng.choice(n_available, size=n_cells, replace=replace)

    X = control_adata.X
    sampled = X[idx].toarray() if hasattr(X, "toarray") else np.asarray(X[idx])
    return sampled


def build_baseline_submission(
    controls_dir: Path,
    targets_csv: Path,
    out_path: Path,
    n_cells: int,
    seed: int,
    cell_type_col: str = "cell_type",
    pert_col: str = "perturbation",
) -> Path:
    pools = load_control_pools(controls_dir)
    targets = pd.read_csv(targets_csv)

    missing_types = set(targets[cell_type_col]) - set(pools.keys())
    if missing_types:
        raise ValueError(
            f"targets_csv references cell types with no control file: {missing_types}"
        )

    rng = np.random.default_rng(seed)
    gene_names = next(iter(pools.values())).var_names.tolist()

    all_X, all_cell_type, all_pert = [], [], []

    for _, row in targets.iterrows():
        cell_type = row[cell_type_col]
        pert = row[pert_col]
        control_adata = pools[cell_type]

        # Sanity check: gene order must match across all control pools,
        # or downstream scoring will silently misalign genes.
        if control_adata.var_names.tolist() != gene_names:
            raise ValueError(
                f"Gene order mismatch in control pool for '{cell_type}' -- "
                "align all control files to the same gene panel/order first."
            )

        pop = sample_baseline_population(control_adata, n_cells, rng)
        all_X.append(pop)
        all_cell_type.extend([cell_type] * n_cells)
        all_pert.extend([pert] * n_cells)

    X = np.vstack(all_X)
    X = np.rint(X).astype(int)  # enforce raw integer counts

    obs = pd.DataFrame({cell_type_col: all_cell_type, pert_col: all_pert})
    adata_out = sc.AnnData(X=X, obs=obs, var=pd.DataFrame(index=gene_names))

    # Spec sanity checks before writing -- catch mistakes here, not at
    # submission time.
    assert (X >= 0).all(), "counts must be non-negative"
    assert X.shape[0] == n_cells * len(targets), "unexpected total cell count"
    expected_per_group = targets.groupby([cell_type_col, pert_col]).size()
    assert (expected_per_group == 1).all(), (
        "targets_csv has duplicate (cell_type, perturbation) rows -- "
        "each pair should appear exactly once."
    )

    adata_out.write_h5ad(out_path)
    print(
        f"Wrote baseline submission: {X.shape[0]} cells x {X.shape[1]} genes "
        f"({len(targets)} (cell_type, perturbation) pairs x {n_cells} cells) -> {out_path}"
    )
    return out_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls-dir", type=Path, required=True)
    parser.add_argument("--targets-csv", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("submission_baseline.h5ad"))
    parser.add_argument("--n-cells", type=int, default=400)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    build_baseline_submission(
        controls_dir=args.controls_dir,
        targets_csv=args.targets_csv,
        out_path=args.out,
        n_cells=args.n_cells,
        seed=args.seed,
    )
