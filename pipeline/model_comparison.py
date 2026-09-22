"""
Model Comparison Harness -- Baseline vs GEARS vs STATE vs Stack
==================================================================

Generalizes validate_on_held_out_h1 (from vcc2026_pipeline.py) to run
several models against the SAME held-out H1 genes and produce one
side-by-side table, so you can see each model's individual
contribution before trusting any ensemble weighting between them.

Kept intentionally simple per your call: this uses the same proxy
metrics as before (Pearson r and sign/direction agreement between
predicted and true log fold change), NOT the real six-metric VCC
scorer. Treat results as a RELATIVE ranking signal between models on
this held-out set, not an absolute estimate of leaderboard score --
the real scorer's perturbation-discrimination and significance-
overlap components in particular don't reduce cleanly to per-gene
correlation.

Each model is wrapped into the same interface:

    predict_fn(train_adata, held_out_genes) -> {gene: lfc_array}

so the harness doesn't need to know anything model-specific.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr

from vcc2026_pipeline import (
    PipelineConfig,
    CONFIG,
    train_gears_on_contexts,
    gears_predict,
    run_state_inference,
)


# ---------------------------------------------------------------------------
# Per-model predict_fn wrappers -- each returns {gene: lfc_array}
# ---------------------------------------------------------------------------

def make_baseline_predict_fn():
    """
    Predicts zero change for every held-out gene -- the same "predict
    no change" floor as the random-control-resampling submission.
    Every other model should be checked against this before you trust
    it's adding real signal.
    """
    def predict_fn(train_adata: sc.AnnData, held_out_genes: list[str]) -> dict[str, np.ndarray]:
        n_genes = train_adata.n_vars
        return {gene: np.zeros(n_genes) for gene in held_out_genes}
    return predict_fn


def make_gears_predict_fn(cfg: PipelineConfig, context_key: str = "cell_type"):
    """
    Trains GEARS once on the training split, reuses it for every held-
    out gene's prediction. Since GEARS wasn't built for cross-cell-type
    transfer, this wrapper assumes train_adata is already restricted to
    contexts appropriate for the target (do that filtering before
    calling the harness, e.g. via rank_donor_contexts).
    """
    def predict_fn(train_adata: sc.AnnData, held_out_genes: list[str]) -> dict[str, np.ndarray]:
        contexts = train_adata.obs[context_key].unique().tolist()
        model = train_gears_on_contexts(train_adata, context_key, contexts)
        return gears_predict(model, held_out_genes)
    return predict_fn


def make_state_predict_fn(model_dir: str, checkpoint: str, control_adata_path: str,
                           pert_col: str = "target_gene", tmp_output: str = "state_tmp_preds.h5ad"):
    """
    Runs STATE's ST module via inference only (pretrained checkpoint,
    no fine-tuning) on the training split's control cells, one call
    covering all held-out genes at once since state tx infer takes a
    perturbation column rather than one gene at a time.
    """
    def predict_fn(train_adata: sc.AnnData, held_out_genes: list[str]) -> dict[str, np.ndarray]:
        preds_adata = run_state_inference(
            model_dir=model_dir,
            checkpoint=checkpoint,
            target_adata_path=control_adata_path,
            pert_col=pert_col,
            output_path=tmp_output,
        )
        control_mean = np.asarray(
            train_adata[train_adata.obs["condition"] == "control"].X.mean(axis=0)
        ).ravel()

        results = {}
        for gene in held_out_genes:
            gene_mask = preds_adata.obs[pert_col] == gene
            pred_mean = np.asarray(preds_adata[gene_mask].X.mean(axis=0)).ravel()
            results[gene] = np.log2((pred_mean + 1e-6) / (control_mean + 1e-6))
        return results
    return predict_fn


def make_stack_predict_fn(model_dir: str, control_adata_path: str,
                           pert_col: str = "target_gene", tmp_output: str = "stack_tmp_preds.h5ad"):
    """
    Stack's whole design point is in-context prediction with no
    fine-tuning: feed it the training split's control cells as context
    and ask for held-out genes directly. Fill in the real Stack CLI/
    Python call here once you've confirmed its actual interface --
    this mirrors run_state_inference's shape since both are inference-
    only forward passes, but Stack's exact call signature hasn't been
    verified against its repo yet.
    """
    def predict_fn(train_adata: sc.AnnData, held_out_genes: list[str]) -> dict[str, np.ndarray]:
        import subprocess
        subprocess.run(
            [
                "stack", "infer",  # TODO: confirm actual Stack CLI subcommand name
                "--model-dir", model_dir,
                "--adata", control_adata_path,
                "--pert-col", pert_col,
                "--targets", ",".join(held_out_genes),
                "--output", tmp_output,
            ],
            check=True,
        )
        preds_adata = sc.read_h5ad(tmp_output)
        control_mean = np.asarray(
            train_adata[train_adata.obs["condition"] == "control"].X.mean(axis=0)
        ).ravel()

        results = {}
        for gene in held_out_genes:
            gene_mask = preds_adata.obs[pert_col] == gene
            pred_mean = np.asarray(preds_adata[gene_mask].X.mean(axis=0)).ravel()
            results[gene] = np.log2((pred_mean + 1e-6) / (control_mean + 1e-6))
        return results
    return predict_fn


# ---------------------------------------------------------------------------
# Comparison runner
# ---------------------------------------------------------------------------

def compare_models(
    h1_adata: sc.AnnData,
    held_out_genes: list[str],
    condition_key: str,
    named_predict_fns: dict[str, callable],
) -> pd.DataFrame:
    """
    Runs each named predict_fn against the same held-out H1 genes and
    returns one long-format DataFrame: one row per (model, gene) with
    pearson_r and direction_agreement, so you can compare distributions
    across models (not just a single mean number that could hide a
    model doing great on some genes and terribly on others).
    """
    train_mask = ~h1_adata.obs[condition_key].isin(held_out_genes)
    train_adata = h1_adata[train_mask]

    control_mean = np.asarray(
        h1_adata[h1_adata.obs[condition_key] == "control"].X.mean(axis=0)
    ).ravel()

    all_rows = []
    for model_name, predict_fn in named_predict_fns.items():
        predictions = predict_fn(train_adata, held_out_genes)

        for gene in held_out_genes:
            true_mean = np.asarray(
                h1_adata[h1_adata.obs[condition_key] == gene].X.mean(axis=0)
            ).ravel()
            true_lfc = np.log2((true_mean + 1e-6) / (control_mean + 1e-6))
            pred_lfc = predictions[gene]

            r, _ = pearsonr(true_lfc, pred_lfc)
            sign_agree = np.mean(np.sign(true_lfc) == np.sign(pred_lfc))
            all_rows.append({
                "model": model_name,
                "gene": gene,
                "pearson_r": r,
                "direction_agreement": sign_agree,
            })

    return pd.DataFrame(all_rows)


def summarize_comparison(results_df: pd.DataFrame) -> pd.DataFrame:
    """
    Per-model summary: mean/median/std of pearson_r and direction
    agreement across held-out genes, plus win-rate against the
    baseline model per gene (fraction of genes where this model beat
    baseline's pearson_r) -- a more robust signal than comparing raw
    means alone when a few genes are outliers.
    """
    summary = results_df.groupby("model").agg(
        pearson_r_mean=("pearson_r", "mean"),
        pearson_r_median=("pearson_r", "median"),
        pearson_r_std=("pearson_r", "std"),
        direction_agreement_mean=("direction_agreement", "mean"),
    ).reset_index()

    if "baseline" in results_df["model"].unique():
        baseline_by_gene = results_df[results_df["model"] == "baseline"].set_index("gene")["pearson_r"]
        win_rates = {}
        for model_name in results_df["model"].unique():
            if model_name == "baseline":
                continue
            model_by_gene = results_df[results_df["model"] == model_name].set_index("gene")["pearson_r"]
            common_genes = model_by_gene.index.intersection(baseline_by_gene.index)
            win_rates[model_name] = (model_by_gene[common_genes] > baseline_by_gene[common_genes]).mean()
        win_rates["baseline"] = np.nan
        summary["win_rate_vs_baseline"] = summary["model"].map(win_rates)

    return summary.sort_values("pearson_r_mean", ascending=False)


# ---------------------------------------------------------------------------
# Example usage sketch
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    cfg = CONFIG

    # h1_adata = sc.read_h5ad("h1_full.h5ad")
    # held_out_genes = ["GENE_A", "GENE_B", "GENE_C"]  # ~20-30 for a meaningful comparison

    # named_fns = {
    #     "baseline": make_baseline_predict_fn(),
    #     "gears": make_gears_predict_fn(cfg),
    #     "state": make_state_predict_fn(
    #         model_dir="/path/to/state_run",
    #         checkpoint="/path/to/state_run/checkpoints/final.ckpt",
    #         control_adata_path="/path/to/h1_controls.h5ad",
    #     ),
    #     "stack": make_stack_predict_fn(
    #         model_dir="/path/to/stack_model",
    #         control_adata_path="/path/to/h1_controls.h5ad",
    #     ),
    # }

    # results_df = compare_models(h1_adata, held_out_genes, "condition", named_fns)
    # summary = summarize_comparison(results_df)
    # print(summary)
    # results_df.to_csv("model_comparison_per_gene.csv", index=False)
    # summary.to_csv("model_comparison_summary.csv", index=False)

    print("Scaffold loaded. Fill in TODOs and the __main__ block with real paths/data.")
