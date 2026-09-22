"""
Virtual Cell Challenge 2026 — Pipeline Framework
==================================================

Synthesizes the approach discussed:
  1. STATE SE embeddings -> identify which donor contexts are closest
     to each unseen target cell line (bridges GEARS's lack of
     cross-cell-type transfer).
  2. Donor effect tables (log fold change) rebuilt from RAW COUNTS
     through one consistent pipeline across public Perturb-seq
     screens (Replogle K562/RPE1, X-Atlas/Orion HCT116/HEK293T,
     Zhu CD4 T-cell resting/stimulated) -- avoids silently mixing
     each paper's own normalization/pseudocount choices.
  3. GEARS trained restricted to (or oversampling) SE-nearest donor
     contexts, since GEARS is not designed for cross-cell-type
     transfer on its own.
  4. STATE ST model run natively on the target line (its designed
     strength: cross-context generalization).
  5. Ensemble of (3) and (4), weighted by validated per-context
     reliability rather than a fixed 50/50 split.
  6. Output as a POPULATION of simulated cells with real variance,
     not a single collapsed point estimate -- a zero-variance
     submission scores far below predicting nothing, since several
     scoring metrics run differential-expression tests that need
     cell-to-cell spread.
  7. A validation harness that holds out real H1 perturbations so
     every component above can be checked against ground truth
     before being trusted on the six actual unseen lines.

Caveats worth remembering while extending this file (from the
conversation, not just implementation detail):
  - Effect transfer across cell types collapses roughly 10x versus
    transfer within the same cell type in different states. Do not
    expect donor-transfer components to carry much signal; weight
    them modestly and validate that assumption on held-out H1
    before trusting it on the real targets.
  - Pure control-cell covariance / linear-response baselines have
    been shown to perform indistinguishably from "predict no
    change" because of signal dilution -- any one perturbation's
    signature is one contributor among many to an aggregate
    covariance. This is why the pipeline leans on GEARS/STATE
    (explicit perturbation conditioning) rather than an unsupervised
    covariance shortcut.

This is a scaffold: it defines the interfaces, data flow, and the
non-obvious correctness constraints. Fill in the marked TODOs with
real calls to `cell-gears`, `arc-state` (the `state` CLI / library),
and your own I/O for each dataset.

CPU-only notes (40 cores, high RAM, no GPU):
  - GEARS training, donor effect tables, and population simulation
    are all CPU-native already -- no change needed, just parallelize
    across cores where the work is embarrassingly parallel (see
    `build_effect_table_parallel` and `simulate_submission_parallel`).
  - SE embedding and Stack/STATE inference are transformer forward
    passes -- CPU-feasible but slower than GPU. Thread env vars below
    and optional dynamic quantization help; batch size should be
    tuned empirically on your hardware, not assumed.
  - Do NOT attempt to train/fine-tune STATE's ST module from scratch
    on CPU -- one reported attempt crashed even on a free GPU Colab
    instance. Use a pretrained checkpoint for inference only, or use
    Stack (in-context, no fine-tuning needed) for the cross-context
    piece instead.
"""

from __future__ import annotations

import os

# Must be set before importing torch/numpy-backed libraries that read
# these at import time. Adjust the "40" if you want to leave headroom
# for the OS / other processes.
_N_CORES = "40"
os.environ.setdefault("OMP_NUM_THREADS", _N_CORES)
os.environ.setdefault("MKL_NUM_THREADS", _N_CORES)
os.environ.setdefault("OPENBLAS_NUM_THREADS", _N_CORES)
os.environ.setdefault("NUMEXPR_NUM_THREADS", _N_CORES)

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from dataclasses import dataclass, field
from pathlib import Path
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr
from concurrent.futures import ProcessPoolExecutor, as_completed

torch.set_num_threads(int(_N_CORES))


# ---------------------------------------------------------------------------
# 0. Config
# ---------------------------------------------------------------------------

@dataclass
class PipelineConfig:
    n_genes: int = 18_533          # 2026 VCC gene panel size
    n_cells_per_pert: int = 400    # submission population size per (context, perturbation)
    cpm_floor: float = 5.0         # genes below this CPM in a screen's controls are zeroed
    # in the donor effect table -- a fold change against a
    # near-zero denominator is dominated by sampling noise.
    se_model_folder: str = "/path/to/SE-600M"
    se_checkpoint: str = "/path/to/SE-600M/se600m_epoch15.ckpt"
    ensemble_weight_gears: float = 0.3   # starting point only -- tune on held-out H1
    ensemble_weight_state: float = 0.5
    ensemble_weight_donor: float = 0.2
    random_seed: int = 0
    n_workers: int = 40           # CPU parallelism for embarrassingly-parallel steps
    quantize_transformers: bool = True  # dynamic quantization for CPU inference


CONFIG = PipelineConfig()


def maybe_quantize(model: torch.nn.Module, cfg: PipelineConfig) -> torch.nn.Module:
    """
    Dynamic quantization for CPU inference on SE/Stack/STATE checkpoints.
    Meaningful speedup on transformer forward passes with minimal
    accuracy loss -- cheap to try, worth timing on a small batch before
    committing to it for a full run.
    """
    if not cfg.quantize_transformers:
        return model
    model.eval()
    return torch.quantization.quantize_dynamic(
        model, {torch.nn.Linear}, dtype=torch.qint8
    )


# ---------------------------------------------------------------------------
# 1. SE embeddings -> context matching
# ---------------------------------------------------------------------------

def compute_se_embeddings(adata_path: str, out_path: str, cfg: PipelineConfig) -> str:
    """
    Wraps `state emb transform`. Run once on a combined h5ad of control
    cells from every context you have (donor screens + H1 + the six
    unseen target lines' controls).

    CLI equivalent:
        state emb transform \
          --model-folder {cfg.se_model_folder} \
          --checkpoint {cfg.se_checkpoint} \
          --input {adata_path} \
          --output {out_path}
    """
    import subprocess
    subprocess.run(
        [
            "state", "emb", "transform",
            "--model-folder", cfg.se_model_folder,
            "--checkpoint", cfg.se_checkpoint,
            "--input", adata_path,
            "--output", out_path,
        ],
        check=True,
    )
    return out_path


def context_vectors(se_adata_path: str, context_key: str = "cell_type") -> pd.Series:
    """Average SE embedding per named context (cell line)."""
    adata = sc.read_h5ad(se_adata_path)
    emb_key = "X_se" if "X_se" in adata.obsm else list(adata.obsm.keys())[0]
    df = pd.DataFrame(adata.obsm[emb_key], index=adata.obs_names)
    df[context_key] = adata.obs[context_key].values
    return df.groupby(context_key).mean().apply(lambda r: r.values, axis=1)


def rank_donor_contexts(
    target_vec: np.ndarray,
    donor_context_vecs: pd.Series,
    top_k: int = 3,
) -> list[tuple[str, float]]:
    """SE-nearest donor contexts for a given (unseen) target line."""
    names = list(donor_context_vecs.index)
    matrix = np.stack(donor_context_vecs.values)
    dists = cdist(target_vec.reshape(1, -1), matrix)[0]
    ranked = sorted(zip(names, dists), key=lambda x: x[1])
    return ranked[:top_k]


# ---------------------------------------------------------------------------
# 2. Donor effect tables (log fold change), rebuilt consistently
# ---------------------------------------------------------------------------

def build_effect_table(
    adata: sc.AnnData,
    condition_key: str,
    control_label: str,
    cfg: PipelineConfig,
    epsilon: float = 1e-6,
) -> pd.DataFrame:
    """
    One consistent LFC estimator applied to every donor screen's RAW
    COUNTS, rather than trusting each paper's own published fold
    changes (different papers use different normalizations and
    pseudocounts -- mixing them silently mixes estimators).

    Returns a (perturbations x genes) DataFrame of log2 fold changes,
    with genes whose control-cell CPM falls below cfg.cpm_floor in
    this screen zeroed out (a ratio against a near-zero denominator
    is noise-dominated, and the challenge's DE-based metrics discard
    those genes anyway).
    """
    counts = adata.X
    total_per_cell = np.asarray(counts.sum(axis=1)).ravel()
    cpm = counts.multiply(1e6 / total_per_cell[:, None]) if hasattr(counts, "multiply") \
        else counts * (1e6 / total_per_cell[:, None])

    control_mask = adata.obs[condition_key] == control_label
    control_mean_cpm = np.asarray(cpm[control_mask.values].mean(axis=0)).ravel()
    low_expr_mask = control_mean_cpm < cfg.cpm_floor

    conditions = [c for c in adata.obs[condition_key].unique() if c != control_label]
    rows = {}
    control_mean_counts = np.asarray(counts[control_mask.values].mean(axis=0)).ravel()

    for cond in conditions:
        pert_mask = (adata.obs[condition_key] == cond).values
        pert_mean_counts = np.asarray(counts[pert_mask].mean(axis=0)).ravel()
        lfc = np.log2((pert_mean_counts + epsilon) / (control_mean_counts + epsilon))
        lfc[low_expr_mask] = 0.0
        rows[cond] = lfc

    return pd.DataFrame(rows, index=adata.var_names).T  # perturbations x genes


def build_effect_table_parallel(
    donor_adatas: dict[str, sc.AnnData],
    condition_key: str,
    control_label: str,
    cfg: PipelineConfig,
) -> dict[str, pd.DataFrame]:
    """
    Runs build_effect_table across donor screens using multiple cores.
    Each screen's table is independent -- ideal fit for your 40 cores.
    AnnData objects are passed by path in real use if they're large
    (pickling a big AnnData across process boundaries is wasteful);
    shown here with in-memory objects for clarity.
    """
    results = {}
    with ProcessPoolExecutor(max_workers=cfg.n_workers) as pool:
        futures = {
            pool.submit(build_effect_table, adata, condition_key, control_label, cfg): name
            for name, adata in donor_adatas.items()
        }
        for fut in as_completed(futures):
            name = futures[fut]
            results[name] = fut.result()
    return results


def combine_effect_tables(
    tables: dict[str, pd.DataFrame],
    weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    """
    Weighted average of per-screen effect tables over their shared
    genes and perturbations. Default weight is uniform; pass in
    validation-derived weights (e.g. inverse SE distance to the
    target context, or reliability measured on held-out H1) instead
    of guessing.
    """
    common_genes = set.intersection(*(set(t.columns) for t in tables.values()))
    common_perts = set.intersection(*(set(t.index) for t in tables.values()))
    weights = weights or {k: 1.0 for k in tables}
    total_w = sum(weights.values())

    combined = None
    for name, table in tables.items():
        sub = table.loc[sorted(common_perts), sorted(common_genes)]
        w = weights[name] / total_w
        combined = sub * w if combined is None else combined + sub * w
    return combined


# ---------------------------------------------------------------------------
# 3. GEARS restricted to SE-nearest contexts
# ---------------------------------------------------------------------------

def train_gears_on_contexts(
    adata_all: sc.AnnData,
    context_key: str,
    contexts_to_use: list[str],
    device: str = "cuda",
):
    """
    Train GEARS restricted to SE-nearest donor contexts. GEARS is not
    designed for cross-cell-type transfer on its own, so this keeps
    it inside contexts an SE-distance check says are plausibly
    related to the target line, rather than asking it to extrapolate
    across all of them at once.
    """
    from gears import PertData, GEARS  # TODO: pip install cell-gears

    subset = adata_all[adata_all.obs[context_key].isin(contexts_to_use)].copy()

    pert_data = PertData("./data")
    pert_data.new_data_process(dataset_name="nn_subset", adata=subset)
    pert_data.load(data_path="./data/nn_subset")
    pert_data.prepare_split(split="simulation", seed=CONFIG.random_seed)
    pert_data.get_dataloader(batch_size=32, test_batch_size=128)

    model = GEARS(pert_data, device=device)
    model.model_initialize(hidden_size=64)
    model.train(epochs=20)
    return model


def gears_predict(model, target_genes: list[str]) -> dict[str, np.ndarray]:
    """Per-gene predicted expression profile, keyed by knocked-down gene."""
    preds = model.predict([[g] for g in target_genes])
    return preds


# ---------------------------------------------------------------------------
# 4. STATE ST inference (native cross-context strength)
# ---------------------------------------------------------------------------

def run_state_inference(
    model_dir: str,
    checkpoint: str,
    target_adata_path: str,
    pert_col: str,
    output_path: str,
    embed_key: str = "X_hvg",
):
    """
    CLI equivalent:
        state tx infer \
          --model-dir {model_dir} \
          --checkpoint {checkpoint} \
          --adata {target_adata_path} \
          --pert-col {pert_col} \
          --embed-key {embed_key} \
          --output {output_path}
    """
    import subprocess
    subprocess.run(
        [
            "state", "tx", "infer",
            "--model-dir", model_dir,
            "--checkpoint", checkpoint,
            "--adata", target_adata_path,
            "--pert-col", pert_col,
            "--embed-key", embed_key,
            "--output", output_path,
        ],
        check=True,
    )
    return sc.read_h5ad(output_path)


# ---------------------------------------------------------------------------
# 5. Ensemble
# ---------------------------------------------------------------------------

def ensemble_lfc(
    gears_lfc: np.ndarray,
    state_lfc: np.ndarray,
    donor_lfc: np.ndarray,
    cfg: PipelineConfig,
) -> np.ndarray:
    """
    Combine the three log-fold-change sources for one perturbation.
    Weights default to cfg values but should be re-fit on the
    held-out-H1 harness below before use -- do not assume 0.3/0.5/0.2
    is correct for your data.
    """
    return (
        cfg.ensemble_weight_gears * gears_lfc
        + cfg.ensemble_weight_state * state_lfc
        + cfg.ensemble_weight_donor * donor_lfc
    )


# ---------------------------------------------------------------------------
# 6. Population generation -- predict a DISTRIBUTION, not one point
# ---------------------------------------------------------------------------

def simulate_perturbed_population(
    control_cells: np.ndarray,       # (n_control_cells, n_genes) raw counts
    predicted_lfc: np.ndarray,       # (n_genes,) log2 fold change vs control mean
    n_cells: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Applies the predicted change to real control cells drawn from the
    target context, rather than fabricating an absolute profile from
    scratch: the control population already encodes which genes are
    on, at what level, and how much cells vary -- for free, and
    without risk of getting that part wrong.

    Critical constraint: DO NOT emit n_cells identical copies of one
    mean profile. Several scoring components run differential-
    expression tests that require population variance; a zero-variance
    submission scores far below even a "predict no change" baseline.
    This samples a *different* subset of control cells per output
    cell, each individually rescaled, to preserve realistic spread.
    """
    n_control = control_cells.shape[0]
    fold = 2.0 ** predicted_lfc  # (n_genes,)

    idx = rng.integers(0, n_control, size=n_cells)
    base = control_cells[idx].astype(np.float64)
    scaled = base * fold[None, :]

    # Binomial resampling to keep realistic per-cell count noise
    # rather than emitting smooth, unrealistically clean floats.
    depth = scaled.sum(axis=1, keepdims=True)
    probs = np.clip(scaled / np.clip(depth, 1e-9, None), 0, 1)
    out = rng.binomial(n=np.clip(depth, 0, None).astype(int), p=probs)
    return np.rint(out).astype(int)


def _simulate_one(args):
    context, gene, controls, lfc, n_cells, seed = args
    rng = np.random.default_rng(seed)
    pop = simulate_perturbed_population(controls, lfc, n_cells, rng)
    return context, gene, pop


def simulate_submission_parallel(
    predictions: dict[str, dict[str, np.ndarray]],
    control_cells_by_context: dict[str, np.ndarray],
    cfg: PipelineConfig,
) -> dict[tuple[str, str], np.ndarray]:
    """
    Parallel version of the per-(context, perturbation) population
    simulation used inside assemble_submission. With six contexts x
    up to a few hundred perturbations each, this is a large but
    fully independent job list -- a good match for 40 cores. Each
    job gets its own seed derived from the base seed so results stay
    reproducible regardless of worker scheduling order.
    """
    jobs = []
    seed_counter = 0
    for context, gene_preds in predictions.items():
        controls = control_cells_by_context[context]
        for gene, lfc in gene_preds.items():
            jobs.append((
                context, gene, controls, lfc,
                cfg.n_cells_per_pert, cfg.random_seed + seed_counter,
            ))
            seed_counter += 1

    results = {}
    with ProcessPoolExecutor(max_workers=cfg.n_workers) as pool:
        for context, gene, pop in pool.map(_simulate_one, jobs, chunksize=4):
            results[(context, gene)] = pop
    return results


# ---------------------------------------------------------------------------
# 7. Validation harness -- hold out real H1 perturbations
# ---------------------------------------------------------------------------

def validate_on_held_out_h1(
    h1_adata: sc.AnnData,
    held_out_genes: list[str],
    condition_key: str,
    predict_fn,  # callable: (train_adata, held_out_genes) -> {gene: predicted_lfc}
) -> pd.DataFrame:
    """
    Cheap, ground-truth-backed sanity check before trusting any
    component (GEARS, donor tables, ensemble weights) on a target
    line where you cannot verify anything. Reports Pearson r between
    predicted and true log fold change per held-out gene, plus
    direction agreement (sign match), which is closer to what the
    real scorer rewards than raw correlation is.
    """
    train_mask = ~h1_adata.obs[condition_key].isin(held_out_genes)
    train_adata = h1_adata[train_mask]

    predictions = predict_fn(train_adata, held_out_genes)

    control_mean = np.asarray(
        h1_adata[h1_adata.obs[condition_key] == "control"].X.mean(axis=0)
    ).ravel()

    rows = []
    for gene in held_out_genes:
        true_mean = np.asarray(
            h1_adata[h1_adata.obs[condition_key] == gene].X.mean(axis=0)
        ).ravel()
        true_lfc = np.log2((true_mean + 1e-6) / (control_mean + 1e-6))
        pred_lfc = predictions[gene]

        r, _ = pearsonr(true_lfc, pred_lfc)
        sign_agree = np.mean(np.sign(true_lfc) == np.sign(pred_lfc))
        rows.append({"gene": gene, "pearson_r": r, "direction_agreement": sign_agree})

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 8. Submission assembly + spec validation
# ---------------------------------------------------------------------------

def assemble_submission(
    predictions: dict[str, dict[str, np.ndarray]],  # {context: {gene: predicted_lfc}}
    control_cells_by_context: dict[str, np.ndarray],
    gene_names: list[str],
    cfg: PipelineConfig,
    out_path: str,
) -> str:
    """
    Builds the final AnnData: raw integer counts, shape
    (n_contexts * n_perturbations * cfg.n_cells_per_pert, cfg.n_genes).
    """
    rng = np.random.default_rng(cfg.random_seed)
    all_X, all_context, all_gene = [], [], []

    for context, gene_preds in predictions.items():
        controls = control_cells_by_context[context]
        for gene, lfc in gene_preds.items():
            pop = simulate_perturbed_population(controls, lfc, cfg.n_cells_per_pert, rng)
            all_X.append(pop)
            all_context.extend([context] * cfg.n_cells_per_pert)
            all_gene.extend([gene] * cfg.n_cells_per_pert)

    X = np.vstack(all_X)
    # Column names/values confirmed against the real `vcc prep` CLI:
    # 'context' (not 'cell_type' or 'context' guessed loosely) and
    # 'target_gene'. Critically, the VALUES in the context column must
    # be exactly the labels used in the downloaded control files (e.g.
    # 'A'/'B'/'C') -- not a filename-derived string like 'context_A'.
    # `all_context` here should already carry those exact labels if it
    # was built from control_cells_by_context keyed the same way
    # baseline_random_controls.py's load_control_pools reads them.
    obs = pd.DataFrame({"context": all_context, "target_gene": all_gene})
    adata_out = sc.AnnData(X=X, obs=obs, var=pd.DataFrame(index=gene_names))
    adata_out.write_h5ad(out_path)

    assert X.shape[1] == cfg.n_genes, f"gene dim mismatch: {X.shape[1]} != {cfg.n_genes}"
    assert X.dtype.kind in "iu" or np.allclose(X, np.rint(X)), "counts must be integers"
    assert (X >= 0).all(), "counts must be non-negative"

    print(f"OK: {X.shape[0]} cells x {X.shape[1]} genes written to {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Example end-to-end sketch (fill in real paths/data before running)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    cfg = CONFIG

    # 1) SE embeddings across every context you have.
    # compute_se_embeddings("all_controls.h5ad", "se_embeddings.h5ad", cfg)
    # ctx_vecs = context_vectors("se_embeddings.h5ad")

    # 2) Rank donor contexts per unseen target line.
    # nearest = rank_donor_contexts(ctx_vecs["target_line_A"], ctx_vecs.drop("target_line_A"))

    # 3) Build + combine donor effect tables (Replogle, X-Atlas/Orion, Zhu CD4T, H1).
    # tables = {name: build_effect_table(adata, "condition", "control", cfg) for name, adata in donor_adatas.items()}
    # donor_table = combine_effect_tables(tables)

    # 4) Train GEARS restricted to SE-nearest contexts.
    # gears_model = train_gears_on_contexts(adata_all, "cell_type", [c for c, _ in nearest])

    # 5) Run STATE ST natively on the target line's controls.
    # state_preds = run_state_inference(...)

    # 6) VALIDATE all of the above on held-out H1 before trusting any of it.
    # report = validate_on_held_out_h1(h1_adata, held_out_genes, "condition", predict_fn)
    # print(report)  # use this to set ensemble weights, not the cfg defaults

    # 7) Ensemble + assemble final submission.
    # assemble_submission(predictions, control_cells_by_context, gene_names, cfg, "submission.h5ad")

    # CPU-specific: prefer running SE / Stack / STATE checkpoints through
    # maybe_quantize(model, cfg) at load time, and prefer the *_parallel
    # variants (build_effect_table_parallel, simulate_submission_parallel)
    # over their single-threaded counterparts wherever the step is
    # embarrassingly parallel across contexts/perturbations.

    print("Scaffold loaded. Fill in TODOs and the __main__ block with real paths/data.")
