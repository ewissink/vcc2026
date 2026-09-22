"""
Baseline submission -- random control resampling
==================================================

For each (context, perturbation) pair required by the challenge,
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
    --controls-dir: one .h5ad per context (e.g. A, B, C), containing
        each containing only that context's control (non-targeting)
        cells, raw counts in .X, genes in .var_names.
    --targets-csv: EITHER
        (a) a flat gene list with a single `target_gene` column --
            the script will cross it with every cell type found in
            controls_dir, since the task is to predict each gene's
            impact in all contexts; OR
        (b) two columns, context,perturbation -- pre-paired rows,
            used as-is.
        Adjust --gene-col / column names below if yours differ.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def load_control_pools(controls_dir: Path, context_col: str = "context") -> dict[str, sc.AnnData]:
    """
    One AnnData of control cells per context label. The label is read
    from each file's OWN obs[context_col] wherever possible -- not
    guessed from the filename -- since the challenge's valid labels
    (e.g. 'A', 'B', 'C') must be reused exactly as given in the
    downloaded control files, and filenames may not match that exactly
    (e.g. a file named context_A.h5ad does NOT mean the label is
    'context_A').

    Falls back to the filename stem, with a 'context_' prefix
    stripped if present, only if the file has no context_col of its
    own -- but this fallback should be treated as a red flag to
    double check against the actual valid label list from `vcc prep`'s
    error output before trusting it.
    """
    pools = {}
    for h5ad_path in sorted(controls_dir.glob("*.h5ad")):
        adata = sc.read_h5ad(h5ad_path)

        if context_col in adata.obs.columns:
            labels = adata.obs[context_col].unique()
            if len(labels) != 1:
                raise ValueError(
                    f"{h5ad_path.name} contains multiple context labels "
                    f"{list(labels)} in obs['{context_col}'] -- expected exactly one "
                    "per control file."
                )
            label = labels[0]
        else:
            label = h5ad_path.stem
            if label.startswith("context_"):
                label = label[len("context_"):]
            print(
                f"Warning: {h5ad_path.name} has no obs['{context_col}'] column -- "
                f"falling back to filename-derived label '{label}'. Verify this "
                "matches the valid labels vcc prep reports before submitting."
            )

        pools[label] = adata
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
    context_col: str = "context",
    pert_col: str = "target_gene",
    gene_col: str = "target_gene",
) -> Path:
    pools = load_control_pools(controls_dir, context_col=context_col)
    targets_raw = pd.read_csv(targets_csv)

    if context_col in targets_raw.columns:
        # Already paired (context, perturbation) rows.
        targets = targets_raw.copy()
        if pert_col not in targets.columns and gene_col in targets.columns:
            targets = targets.rename(columns={gene_col: pert_col})
    else:
        # Flat gene list (e.g. a single `target_gene` column): cross with
        # every context found in controls_dir, since the task is to
        # predict each gene's impact in all contexts, not a pre-paired
        # subset. This is the common VCC shape -- one target-gene list,
        # applied identically across the validation/test cell lines.
        if gene_col not in targets_raw.columns:
            raise ValueError(
                f"targets_csv has neither '{context_col}' nor '{gene_col}' column -- "
                f"got columns: {list(targets_raw.columns)}"
            )
        genes = targets_raw[gene_col].unique()
        contexts = sorted(pools.keys())
        targets = pd.DataFrame(
            [(ct, g) for ct in contexts for g in genes],
            columns=[context_col, pert_col],
        )
        print(
            f"Expanded {len(genes)} genes x {len(contexts)} contexts "
            f"({contexts}) -> {len(targets)} (context, perturbation) pairs"
        )

    missing_contexts = set(targets[context_col]) - set(pools.keys())
    if missing_contexts:
        raise ValueError(
            f"targets reference contexts with no control file: {missing_contexts}. "
            f"Valid contexts from control files: {sorted(pools.keys())}"
        )

    rng = np.random.default_rng(seed)
    gene_names = next(iter(pools.values())).var_names.tolist()

    all_X, all_context, all_pert = [], [], []

    for _, row in targets.iterrows():
        context = row[context_col]
        pert = row[pert_col]
        control_adata = pools[context]

        # Sanity check: gene order must match across all control pools,
        # or downstream scoring will silently misalign genes.
        if control_adata.var_names.tolist() != gene_names:
            raise ValueError(
                f"Gene order mismatch in control pool for '{context}' -- "
                "align all control files to the same gene panel/order first."
            )

        pop = sample_baseline_population(control_adata, n_cells, rng)
        all_X.append(pop)
        all_context.extend([context] * n_cells)
        all_pert.extend([pert] * n_cells)

    X = np.vstack(all_X)
    X = np.rint(X).astype(int)  # enforce raw integer counts

    # Column names confirmed against the real `vcc prep` CLI: 'context'
    # (not 'cell_type') and 'target_gene', with context VALUES that must
    # exactly match the labels the control files themselves use (e.g.
    # 'A'/'B'/'C') -- not filename-derived strings.
    obs = pd.DataFrame({context_col: all_context, pert_col: all_pert})
    adata_out = sc.AnnData(X=X, obs=obs, var=pd.DataFrame(index=gene_names))

    # Spec sanity checks before writing -- catch mistakes here, not at
    # submission time.
    assert (X >= 0).all(), "counts must be non-negative"
    assert X.shape[0] == n_cells * len(targets), "unexpected total cell count"
    expected_per_group = targets.groupby([context_col, pert_col]).size()
    assert (expected_per_group == 1).all(), (
        "targets_csv has duplicate (context, perturbation) rows -- "
        "each pair should appear exactly once."
    )

    adata_out.write_h5ad(out_path)
    print(
        f"Wrote baseline submission: {X.shape[0]} cells x {X.shape[1]} genes "
        f"({len(targets)} (context, perturbation) pairs x {n_cells} cells) -> {out_path}"
    )
    return out_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls-dir", type=Path, required=True)
    parser.add_argument("--targets-csv", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("submission_baseline.h5ad"))
    parser.add_argument("--n-cells", type=int, default=400)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--gene-col", type=str, default="target_gene",
                         help="Column name in targets_csv for a flat gene list (ignored if the CSV already has a context column)")
    args = parser.parse_args()

    build_baseline_submission(
        controls_dir=args.controls_dir,
        targets_csv=args.targets_csv,
        out_path=args.out,
        n_cells=args.n_cells,
        seed=args.seed,
        gene_col=args.gene_col,
    )
