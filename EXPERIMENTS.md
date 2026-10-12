# Experiment log

Running record of submissions, local diagnostics and what they imply. Newest decisions at the top
of each section. Leaderboard scores are SCALED (0 = mean-response baseline, 1 = replicate) unless
marked raw.

## Current plan (2026-10-12)

1. Submit v5 / v6 (donor only, PIE switched off) to test whether cross-line transfer from
   K562 GWPS gives PDS on the real panel, as it does on the H1 proxy fold.
2. Rebuild the local harness with `--panel-match` (effect sizes matched to the Challenge panel)
   and run the 09 grid on it once the `pie prep` command for the harness context is written.
3. Fix the shared response (nmae -0.13 in v4) on the harness.
4. Optional: per-pert donor weight in 09 from `split_half_challenge.csv` (job still running).

Open items: v1-v3 scores and flags (not accessible right now); `pie prep` command used for
`vcc2026_val`.

Context identities (expected, from CCLE matching): **A = Jurkat, B = HeLa, C = CAL33 or CAL27**.
No public Perturb-seq screen found for HeLa, CAL27 or CAL33 (web search 2026-10-12). Genome-wide
donors that cover the panel: K562 GWPS (have), X-Atlas/Orion HCT116 and HEK293T (not fetched),
KOLF2.1J iPSC CRISPRi (figshare 10.25452/figshare.plus.27261219.v1, not fetched).

## Submissions

All on partition `val`, panel `vcc2026-val-1`, anchors `vcc2026-valA-r4+valB-r4+valC-r4`.

| Version | 09 flags (beyond --pred/--pie-meta/--data-dir/--pert-counts) | Overall | PDS | MSE | JAC | NMAE | FID | REACH | Rank |
|---|---|---|---|---|---|---|---|---|---|
| v1-v3 | TODO | | | | | | | | |
| v4 `pred_v4_k562_0.1` | `--shrink 0 --donor-npz replogle_k562_gwps.npz --donor-scale 0.1` | 0.0106 | 0.030 | 0 | 0.010 | -0.133 | -0.028 | 0.184 | 953 |
| v5 `pred_v5_donly_0.1_nc` | `--shrink 0 --lfc-scale 0 --donor-npz replogle_k562_gwps.npz --donor-scale 0.1 --no-donor-center` | -0.242 | 0.125 | 0 | -0.070 | 0.022 | **-1.562** | 0.037 | 1292 |
| ~~v6 `pred_v6_donly_0.3_nc`~~ | same, `--donor-scale 0.3` | dropped after v5 (FID collapse) | | | | | | | |
| v6 `pred_c_pie0.5_d0.3` | `--shrink 0 --lfc-scale 0.5 --donor-npz replogle_k562_gwps.npz --donor-scale 0.3` | **0.0267** | 0.102 | 0 | 0.008 | -0.128 | -0.022 | 0.201 | 926 |

v6 vs v4: best so far. PDS 0.03 -> 0.10, REACH up slightly, FID unchanged at 450-1,050 DE calls
per pert (so N_conf is probably below ~500-1,000), but NMAE did not move (-0.128) despite halving
PIE's shared part -- so NMAE is not simply PIE overshooting, or donor 0.3 added back what PIE 0.5
removed (v5, with no PIE and donor 0.1, had NMAE +0.02). PDS from the K562 donor on the real panel
(0.10-0.125 scaled across v5/v6) is far below what the H1 proxy fold suggested.

v6 raw: PDS 0.546, MSE 13.619, JAC 0.033, NMAE 1.071, FID 0.506, REACH 0.256.
- Raw NMAE > 1 = worse than predicting no change on the true DE genes (proxy no_change = 1.0;
  leaderboard rank 900 = 0.979). Raw MSE 13.6 is ~4x rank 900 (3.4). Together with FID ~0.5
  (coin-flip directions), PIE's shared part looks like direction noise: needed only to produce
  DE calls, and it costs NMAE and MSE.
- Raw PDS 0.546: K562 transfer on the real panel is weaker than even the K562 proxy fold (0.59).

v5 vs v4 (same 0.1 x centered per-pert donor part; only the shared part differs: PIE panel mean
in v4, 0.1 x K562 panel mean in v5): smaller shared part -> PDS 0.03 -> 0.125 and NMAE
-0.13 -> +0.02, but FID collapsed (-1.56) and REACH fell 0.18 -> 0.04. So PIE's shared part,
though too large, carried the right directions; the K562 mean response (or having no good shared
part) does not.

v5 diagnostics: (1) generated cells carry the intended directions (sign agreement with the K562
delta on each pert's top-200 genes 0.78-0.83; library-size ratio ~1.00; median shift ~0), so it is
not a cell-generation artifact. (2) K562 panel mean vs PIE panel mean per context: Pearson
-0.02 / 0.01 / 0.01 (A/B/C), top-500 sign agreement ~0.59-0.62 -- unrelated, not opposite. Neither
explains a raw FID far below 0.5 by itself; need v5's raw values.

**Resolved: v5 raw FID = 0.048 is a YIELD collapse, not wrong directions.** The scored member is
`de_wilcoxon_direction_fidelity_yield_raw = k / max(n_pred, N_conf)` (cell-eval2
`metrics/direction.py:776`; base 0.5 = coin-flip sign): calling fewer DE genes than the truth caps
the score at ~n_pred/N_conf. Above N_conf it is roughly sign accuracy over the called genes (v4:
~0.5, i.e. PIE's shared directions are coin-flip quality; top teams only reach ~0.56). REACH is
probably also count/rank based (purity floor 0.9 on genes ranked by significance).

Local DE-call counts (Wilcoxon vs 2,000 context controls, BH q < 0.05, genes > 5 CPM, 25 perts per
context; approximates but is not cell-eval2's DE):

| File | median n_sig per pert | A | B | C |
|---|---|---|---|---|
| v4 (PIE 1.0 + donor 0.1 centered) | 4379 | 4379 | 5049 | 2335 |
| v5 (donor 0.1 uncentered only) | 1 | 1 | 1 | 0 |
| PIE 0.3 + donor 0.1 | 42 | 26 | 68 | 42 |
| PIE 0.5 + donor 0.3 | 739 | 644 | 1056 | 454 |
| PIE 1.0 + donor 0.3 | 4193 | 4193 | 4940 | 2403 |

=> the PIE shared scale drives the number of DE calls, steeply. PIE 0.3 would have collapsed
like v5 (not submitted). Rule before any submission: run the n_sig count; reject anything far
below v4's calls until N_conf is known.

v4 reading: with `--shrink 0` every pert got PIE's panel-mean effect (identical across perts)
plus 0.1 x centered GWPS delta. The shared PIE part dominated each prediction's direction, so
perts were barely distinguishable (PDS ~0.52 raw), and it is worse than the true mean response
(NMAE, FID negative). v5/v6 drop the PIE part (`--lfc-scale 0`) and keep the donor's own
panel-mean response (`--no-donor-center`), matching the proxy's `transfer` rows. Watch: PDS
should rise; REACH/FID may fall (the mean response carries most of those).

Commands:

```bash
python scripts/09_make_cells.py --pred pred_val_full.parquet --pie-meta /workdir/emw97/pie/data/vcc2026_val/meta.json --data-dir data_dir --pert-counts data_dir/pert_counts.csv --out pred_v5_donly_0.1_nc.h5ad --workers 32 --shrink 0 --lfc-scale 0 --donor-npz results/harmonized/replogle_k562_gwps.npz --donor-scale 0.1 --no-donor-center
python scripts/09_make_cells.py --pred pred_val_full.parquet --pie-meta /workdir/emw97/pie/data/vcc2026_val/meta.json --data-dir data_dir --pert-counts data_dir/pert_counts.csv --out pred_v6_donly_0.3_nc.h5ad --workers 32 --shrink 0 --lfc-scale 0 --donor-npz results/harmonized/replogle_k562_gwps.npz --donor-scale 0.3 --no-donor-center
```

## Leaderboard reference points (2026-10-12, scaled with raw in parentheses)

| Rank | Overall | PDS | MSE | JAC | NMAE | FID | REACH |
|---|---|---|---|---|---|---|---|
| 1 | 0.463 | 0.879 (0.896) | 0.776 (0.251) | 0.076 (0.059) | 0.431 (0.752) | 0.174 (0.563) | 0.457 (0.483) |
| 100 | 0.254 | 0.745 (0.836) | 0.329 (0.675) | 0.003 (0.031) | 0.175 (0.894) | 0.045 (0.526) | 0.229 (0.282) |
| 500 | 0.153 | 0.569 (0.769) | 0 (5.093) | -0.10 (0.027) | 0.155 (0.907) | 0.006 (0.514) | 0.170 (0.230) |
| 900 | 0.035 | 0.259 (0.619) | 0 (3.436) | -0.20 (0.023) | 0.033 (0.979) | -0.155 (0.478) | 0.053 (0.126) |

- PDS is the main gap to the top 100 (0.03 vs 0.75 is ~+0.12 overall); NMAE second (~+0.05).
- Scaled MSE is floored at 0 below roughly the top 100. JAC is ~0 for everyone.
- Overall is close to, but not always exactly, the mean of the six scaled metrics.

## Scoring rules (cell-eval2 0.15.0, `rule_version 3`)

From comments in `cell_eval2/competition.py` and `baseline.py`:
- `pds_cosine` ranks with EVERY panel target gene removed (#343). Spiking target genes now scores
  the control-paste floor (raw 0.5000). The mean-response baseline's raw PDS on the official
  contexts is ~0.53 / 0.53 / 0.51.
- Three scored vcc2026 metrics no longer score each pert's own target gene (#172).
- => Explicit target knockdown does not help. PDS needs downstream (off-target) effects.
- Our `04_score_baselines.py` proxy also removes all panel targets, so it matches on this point.

## Diagnostics

### PIE checkpoint `replogle_xdataset` training data
`pie_replogle_nadig_essential`, `pie_tahoe100m`, `pie_jiang`, `pie_arc_vcc_25`,
`pie_x_atlas_orion` (config also has `dataset_weights: replogle: 0.0`, not yet checked).
K562 GWPS is NOT in it, so harness perts (essential-screen perts excluded) are unseen in K562.

### GWPS donor coverage of the Challenge panel
269/300 Challenge perts are in K562 GWPS. Their GWPS delta norm (own target included): median
14.4, IQR 12.4-17.1 -- the same as all GWPS perts (median 14.1), i.e. typical, not enriched.

### Transfer scale alpha (04, leave-one-cell-line-out)
| Fold | alpha_other | alpha_self |
|---|---|---|
| K562 | 0.60 | 0.30 |
| HepG2 | 0.45 | 0.61 |
| Jurkat | 0.48 | 0.52 |
| RPE1 | 0.43 | 0.74 |
| H1 (10x Flex) | 0.49 | 0.07 |

### 04 proxy tables (raw; pds 0.5 = no skill; lfc_nmae lower is better)
| Method | K562 pds | K562 reach | K562 nmae | H1 pds | H1 reach | H1 nmae |
|---|---|---|---|---|---|---|
| oracle_mean | 0.500 | 0.192 | 0.987 | 0.500 | 0.268 | 0.956 |
| transfer | 0.544 | 0.162 | 0.956 | 0.797 | 0.110 | 0.971 |
| transfer_scaled_other | 0.593 | 0.102 | 0.958 | 0.814 | 0.109 | 0.973 |
| transfer_topn | 0.597 | 0.034 | 0.982 | 0.847 | 0.037 | 0.976 |

K562 fold: 5 panels of 300, `stratum all` (novel perts sit at 0.5). H1 fold: one panel of all
210 VCC-2025 H1 perts, mostly covered by GWPS -- the closest analogue to the Challenge (Flex,
distant context, genome-wide donor). Transfer there reaches ~rank-100 raw PDS. Caveats: single
panel, H1 perts possibly effect-enriched, pseudobulk proxy.

### Harness `harness_K562gw` (11, top-N panel)
7,679 / 18,533 Challenge genes in GWPS; 4,864 qualifying perts; panel = top 300 by effect norm
(21.1-31.0 vs median 14.1); 37,664 context + 37,664 held-out controls; ~67 cells per panel pert.
- Spearman(cells, effect norm) = -0.546.
- Split-half cosine (all genes, own target dropped): median 0.044; by cell tertile 0.020 / 0.024 /
  0.109. => the top-N panel is largely selected on noise. Superseded by `--panel-match`:

```bash
python scripts/11_build_harness_data.py --h5ad perturb_datasets/ReplogleWeissman2022_K562_gwps.h5ad --harmonized-npz results/harmonized/replogle_k562_gwps.npz --gene-names data_dir/gene_names.csv --name K562 --out harness_K562gw_match --panel-match data_dir/pert_counts.csv --exclude-npz results/harmonized/replogle_k562.npz results/harmonized/replogle_rpe1.npz results/harmonized/nadig_jurkat.npz results/harmonized/nadig_hepg2.npz
```

### GWPS split-half on Challenge perts (running)
Output `split_half_challenge.csv`: per Challenge pert, cells / effect norm / split-half cosine /
top-200-gene sign agreement, plus 100 control-cell null groups. Key number: how many Challenge
perts beat the null's 95th percentile of split-half cosine.
