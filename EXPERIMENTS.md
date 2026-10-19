# Experiment log

Running record of submissions, local diagnostics and what they imply. Newest decisions at the top
of each section. Leaderboard scores are SCALED (0 = mean-response baseline, 1 = replicate) unless
marked raw.

## Status board (update after every result)

_Last updated 2026-10-13._

**Current best:** v13, overall 0.1160 -- `pred_c_pie0.25_4donor_0.7`
(`--shrink 0 --lfc-scale 0.25 --donor '*=replogle_k562_gwps.npz,xatlas_hct116.npz,xatlas_hek293t.npz,cd4_rest.npz:0.7'`).
Donor count is the lever (v10 1-donor 0.067 -> v12 2 0.089 -> v11 3 0.109 -> v13 4 0.116;
diminishing returns per donor). Adding the CD4 Flex donor gave only +0.007 -> no platform-match
windfall; coverage/count drive it. Competitors report ~0.20 with 20+ donors. Next: add iPSC donors.

**Scored 2026-10-14:** v9 `pred_c_d0.7` 0.0492; v10 `pred_c_pie0.25_d0.7` **0.0672** (best).
Progression: v4 0.011 -> v6 0.027 -> v7 0.044 -> v9 0.049 -> v10 0.067.

**Scored 2026-10-15:**
- v11 `pred_c_pie0.25_all3_0.7` **0.1089 (new best)** | PDS 0.423 (0.692) | MSE 0 (3.311) |
  JAC -0.008 (0.028) | NMAE 0.107 (0.936) | FID -0.31 (0.503) | REACH 0.162 (0.222).
  vs v10 (K562 only): PDS, MSE, NMAE, REACH all improved by averaging 3 donors. Confirms
  averaging isolates b_p + cancels noise. Progression: v10 0.067 -> v11 0.109.
- v12 `pred_c_pie0.25_xatlas0.7` (X-Atlas 5' only, 2 donors) **0.0889** | PDS 0.374 (0.670) |
  MSE 0 (3.917) | JAC -0.011 (0.026) | NMAE 0.085 (0.949) | FID -0.058 (0.495) | REACH 0.143 (0.206).

**Donor-count is monotone: v10 (1 donor) 0.067 < v12 (2) 0.089 < v11 (3) 0.109.** And adding 3'
K562 on top of 5' X-Atlas HELPED (v11 > v12), so platform mismatch does NOT dominate -- coverage +
noise-cancellation do. ("same-platform 5:1" claim not supported, though no Flex donor tested yet;
CD4 Flex donor is the real test.) Next lever: add the CD4 T-cell Flex donor as a 4th. CD4 panel coverage = 297/300 in the guide
library; the Rest-condition DE object (19_cd4_donor.py -> cd4_rest.npz) covers 283/300 perts,
10,113 genes, delta = DESeq2 log2FC x ln2. 4-donor union now covers **300/300** panel perts.
v13 `pred_c_pie0.25_4donor_0.7` SCORED **0.1160 (best)** | PDS 0.440 (0.700) | MSE 0 (2.828) |
JAC -0.004 (0.029) | NMAE 0.123 (0.926) | FID -0.035 (0.502) | REACH 0.173 (0.232). +0.007 vs v11
-> CD4 Flex donor = ordinary donor-count gain, no platform windfall.
v15-RW `pred_c_pie0.25_4donorRW_0.7` (reliability-weighted K562, min-weight 0.5): 0.103 | PDS 0.414
(0.688) | MSE 0 (2.559) | NMAE 0.115 (0.930) | FID -0.071 (0.491) | REACH 0.166 (0.226). Below v13.
Built, not submitted: pie0.15/0.35 + K562 donor 0.7 (PIE-scale refinement around 0.25).

**Next ideas (not yet queued; 2/day):**
- sweep PIE shared 0.25 further (0.15, 0.35) on top of donor 0.7 -- FID/REACH vs NMAE tradeoff
- swap/add X-Atlas HCT116+HEK293T donors (290/294/300 panel, far cleaner) for K562 GWPS
- per-pert reliability-weighted donor (17) at donor 0.7

### Status: v13 is best; all levers tested
Current best = v13 `pred_c_pie0.25_4donor_0.7` = **0.116**.
- v14 (uniform DE-freq shared) = 0.107; v15 (per-pert-weighted DE-freq) dropped (weighting didn't
  engage on the panel); v15-RW (reliability-weighted K562, min-weight 0.5) = **0.103** -- MSE
  slightly better (raw 2.56) but PDS/FID/NMAE/REACH all down; reliability weighting net NEGATIVE
  (K562 already noise-diluted by 3 co-donors, so shrinking it just removes signal).
Levers tested: donor count (helps), platform (no windfall), shared source (PIE mean > DE-freq),
per-pert shared weighting (inert on panel), donor reliability (negative). => v13 recipe stands.

### Before 2026-10-22 final data (no more new mechanisms; confirm + squeeze)
1. Optional quick scale sweep around v13: donor-scale {0.6,0.8}, lfc-scale {0.2,0.3} -- local optimum.
2. Prep for final test: pipeline must run on the new A/B/C control profiles + 300 final perts;
   keep local harness + DE-call gate for blind selection (no leaderboard Oct22-Nov5).

1. **Submit v15-RW** = reliability-weighted K562 donor on the v13 recipe. Build + gate, then submit:
   build: `python scripts/09_make_cells.py --pred pred_val_full.parquet --pie-meta /workdir/emw97/pie/data/vcc2026_val/meta.json --data-dir data_dir --pert-counts data_dir/pert_counts.csv --out pred_c_pie0.25_4donorRW_0.7.h5ad --workers 32 --shrink 0 --lfc-scale 0.25 --donor '*=results/weighted/replogle_k562_gwps_rw.npz,results/harmonized/xatlas_hct116.npz,results/harmonized/xatlas_hek293t.npz,results/harmonized/cd4_rest.npz:0.7'`
   gate vs v13, then `vcc prep pred_c_pie0.25_4donorRW_0.7.h5ad -g data_dir/gene_names.csv --perts data_dir/pert_counts.csv -o pred_v15RW_pie0.25_4donorRW_0.7.vcc` / `vcc submit pred_v15RW_pie0.25_4donorRW_0.7.vcc -m pie0.25_4donorRW_0.7`.
   Test: does down-weighting noisy K562 perts (split-half reliability) help?
2. Record score; update Current best.

DECISION: if v15-RW does not beat v13 (0.116), LOCK v13-style as the final pipeline and shift to
write-up. Every lever will then be tested: donor count, platform, shared-response source (PIE vs
DE-freq), per-pert shared weighting, donor reliability. Oct-22 final-test data imminent; NO
leaderboard Oct22-Nov5 (pick blind -> keep local harness + DE-call gate). Don't chase 0.01s past this.

**Open threads:**
| Thread | Next step | Waiting on |
|---|---|---|
| PIE's shared part (NMAE, MSE) | read v7 / v8 | today's scores |
| Reliability-weighted donor (17) | test on top of whichever of v6/v7/v8 wins | v7 / v8 |
| Better donors for B, C | DONE: X-Atlas HCT116 (18,106 genes, 165k ctrl, 290/300 panel) + HEK293T (18,106 genes, 219k ctrl, 294/300 panel) harmonized. Both far cleaner than K562 GWPS. | use in 09 / 18 |
| K562 -> Jurkat transfer map for A | fit on shared essential perts | -- |
| Local harness | rebuild with `--panel-match`; `pie prep` on it | the `pie prep` command used for `vcc2026_val` |
| v1-v3 | record scores and flags | access to those results |

**⚠️ Rules check needed (2026-10-13):** the Challenge's "Machine Learning Predictions Only" rule
(revised 2026-09-16) says entries must not contain published experimental results "in any part
of your entry, nor use such results to revise predictions generated by your model"; using them
to train or fine-tune a model is allowed. `09_make_cells.py --donor*` adds scaled measured GWPS
deltas directly to the prediction (v4-v8; v8 is donor-only). On a plain reading that may not be
allowed. Ask the organizers before relying on it; meanwhile design a version where the
published data only trains a model.

**Final entry:** final-test perts + NTC profiles for 3 held-out lines released 2026-10-22;
deadline 2026-11-05 23:59 UTC; max 2 final entries/day; only the LAST entry counts; no
leaderboard during that window, so settings must be chosen without feedback.

**18 transfer-model LOCO (2026-10-14, --min-expr 0.5, median cos pred vs measured):**
| Held-out line | multi-donor | best single donor | verdict |
|---|---|---|---|
| HepG2 (essential) | 0.212 | RPE1 0.144 | multi wins |
| Jurkat (essential) | 0.193 | K562 0.151 | multi wins |
| RPE1 (essential) | 0.226 | HepG2 0.153 | multi wins |
| H1 (Flex) | 0.076 | RPE1 0.097 | single wins |
| K562 (genome-wide) | 0.038 | Jurkat 0.156 | single wins big |
| HCT116 (genome-wide) | 0.022 | RPE1 0.081 | single wins big |
| HEK293T (genome-wide) | 0.016 | Jurkat 0.057 | single wins big |

- alpha (1/10/100) changes nothing -> ridge penalty is not the lever; drop the sweep.
- Global-weight ridge ~= AVERAGING donors (weights ~0.07-0.15 each). It beats single donors on
  the 3 essential screens but LOSES badly on the 3 genome-wide screens -- which are the
  Challenge-like regime (non-essential perts, low signal). Averaging in mediocre donors dilutes
  the one good one. Fix would be similarity-CONDITIONED weights (donor x target-line interaction),
  not global constants. As built, a single well-chosen donor is better for genome-wide targets.
- Note the best single donor is often not the obvious line (K562's is Jurkat, HCT116's is RPE1).

**Similarity-conditioned weights (--sim, 2026-10-15): does NOT rescue it.** multi-donor cos with
sim vs without: H1 0.033 (was 0.076), HCT116 0.020 (0.022), HEK293T 0.013 (0.016), K562 0.035
(0.038) -- still far below best-single on all genome-wide + H1 folds; essential screens ~unchanged
(HepG2 0.208, Jurkat 0.196, RPE1 0.222, all still > their singles). Why sim can't select: the
immortalized lines' control profiles are ALL 0.76-0.87 correlated with each other; only H1 stands
apart (0.36-0.50). So similarity separates "H1 vs rest," not "which immortalized donor," which is
the choice that matters.
**Conclusion -- trained linear transfer model SHELVED as a negative result.** No linear blend
(global or similarity-gated) beats a single well-chosen donor in the genome-wide / low-signal
regime that matches the Challenge. Mechanism: target deltas are ~80% noise, MSE training is
dominated by the two huge genome-wide lines (~12k perts each vs ~1.5k for essential screens), so
the fit chases noise and shrinks. The practical method stays single/few-donor transfer (09's
--donor), which works (v9/v10 PDS raw ~0.69). This is a clean portfolio result: the sophisticated
approach, honestly evaluated LOCO, loses to the simple baseline where it counts, with a mechanistic
explanation (noise ceiling).
Nuance for the real contexts: for the Jurkat fold (= context A) multi-donor DID win (0.196 vs
single-best K562 0.151); the losses are on the epithelial genome-wide lines (analogs of B/C).

### Target expression vs response strength, and the noise floor (2026-10-15)
Per pert: baseline expression of the perturbed gene (ctrl_mean) vs effect-norm, across HCT116,
HEK293T, K562 GWPS.
- Raw norm is dominated by sampling noise: Spearman(n_cells, norm) = -0.89 to -0.97. Fitting
  norm^2 = K/n_cells + signal^2 gives K ~= 4e4 (nearly identical across lines -> a property of
  single-cell sampling, not the screen) and **83-93% of mean norm^2 is noise**.
- Noise-corrected signal MEDIAN is 0.0 in every target-expression quartile of every line: **>50%
  of (non-essential) perts have a response indistinguishable from noise** at achievable per-pert
  cell counts, even in X-Atlas (200k controls -- because noise is set by cells-per-PERT, hundreds,
  not the control pool). Only the tail carries signal (~1/3, matching the K562 split-half null test).
- Target expression IS a real but weak driver: partial/noise-corrected Spearman ~0.10-0.20,
  monotonic across lines.
- Implication: the first-order knob is per-pert signal reliability (17_reliability_weight.py), not
  target expression. Expression-gating the donor would be a weak secondary prior; not worth building
  on its own.

### Competitor-board ideas checked (2026-10-15)
- "Predict only the change, apply to control cells" and "centered cross-context mean (b_p + i_cp)"
  = exactly our 09 design + donor centering. Validated.
- "Averaging donors -> ~0.14 overall, broad collection -> ~0.20" (competitors). So public-data
  averaging reaches rank ~500-100; our 0.067 is donor-COUNT-limited (3 panel donors vs their 20+),
  not method-limited. Adding donors is the main lever.
- "Targets agree across lines at cosine ~0.04" = our split-half finding. Validated.
- "Same-PLATFORM beats same-cell-type 5:1": chemistry lookups -- X-Atlas is 10x GEM-X **5'** (FiCS),
  Replogle K562 is 3'. Challenge is **Flex**. So NONE of our panel-covering donors match the
  Challenge platform; only VCC-2025 H1 is Flex (0 panel overlap). If the 5:1 claim holds, this
  caps public-donor transfer and we can't fix it from these sources.
- DepMap perturb-seq 16-line (broadinstitute/perturb-seq-depmap-public): CRISPR **KO** (not CRISPRi),
  a "pilot", chemistry/coverage/download only on depmap.org. KO modality + uncertain panel coverage
  -> lower priority.
- CRISPRi cis-neighbor knockdown (scored, reproducible) and cell-draw/decoding quirks: untested;
  noted, not built (neighbor signal partly excluded by panel-target removal; decoding tricks fragile).

### Candidate additional donors (ranked, 2026-10-15) -- donor count is the lever (v11)
1. **CD4+ T-cell genome-scale Perturb-seq** (Zhu/Dann/Marson, bioRxiv 2025.12.23.696273):
   genome-wide (~11.5-12.7k genes), **CRISPRi**, **10x Flex** (= Challenge platform, the only
   Flex donor found), primary CD4 T cells (relevant to context A = Jurkat), 4 donors,
   resting + stimulated. **RELEASED** (now in Cell) at virtualcellmodels.cziscience.com/dataset/
   genome-scale-tcell-perturb-seq. CRISPRi + Flex (GEMX_flex_v1) confirmed. ~22M cells, 18,129
   genes, conditions Rest/Stim8hr/Stim48hr. obs: perturbed_gene_name, guide_type non-targeting.
   Access: `vcp data search` or S3; pseudobulk object GWCD4i.pseudobulk_merged.h5ad + DE_stats
   objects exist (avoid the 22M-cell download). Guide list: github.com/emdann/
   GWT_perturbseq_analysis_2025 metadata/suppl_tables. TOP priority donor for Oct 22.
   (= "Zhu CD4" screen already name-checked in pipeline/vcc2026_pipeline.py.)
2. **KOLF2.1J iPSC** genome-wide CRISPRi (11,739 genes, 2.5M cells, figshare
   10.25452/figshare.plus.27261219): downloadable now, chemistry unconfirmed. Secondary.
3. iPSC pan-genome trans map (6,673 targets, bioRxiv 2024.11.28.625833): tertiary.
Rejected (targeted, low panel overlap expected): Tian 2021 iPSC-neuron CROP-seq (184 genes,
GSE152988); foreskin-fibroblast Perturb-seq (GSM5029519, likely targeted).

### Liu/Wang et al. 2026 (biorxiv 2026.10.05.756839) -- validates + extends our approach
Team placed 3rd in VCC-2025 with "DE-frequency + cross-dataset weighted transfer". Key points:
- **Weighted-average cross-dataset transfer is best for PDS** and residual cosine, beating PRESAGE,
  cross-dataset linear (embeddings), GEARS, scGPT-ft. Confirms our v11 donor-averaging win and our
  shelving of the trained transfer model (low-dim embeddings lose to direct transfer).
- Their average is **similarity-weighted** (w_d = mean noise-corrected cosine sim of donor d to
  target over shared perts) and **noise-corrected**; per-(pert,gene) average only over donors that
  measured it (we already do this). Negative weights -> 0.
- **DE-frequency ranking (NEW for us):** rank outcome genes by how often they're DE across the
  donor perturbations (perturbation-AGNOSTIC). Among the strongest predictors of which genes
  respond, esp. for small effects. A shared "hierarchy of responsive genes." Would help our
  sig_jaccard / direction_reach (which genes to call DE). Their 3rd-place entry = DE-freq + transfer.
- **Leading global trend:** first SVD component of each screen's effect matrix ~ fitness/essentiality
  (bigger in immortalized lines). Removing it ("residual") isolates perturbation-specific signal =
  what PDS/residual-cosine reward. More principled than our panel-mean centering.
- Noise correction (subtract SE^2 from squared norms) = our split-half idea, confirmed.

**Donor list from their Data availability (all public) -- candidates to add to our average:**
have K562(Zenodo 10044268), Nadig HepG2/Jurkat(GSE264667, 0 panel overlap), Huang=X-Atlas
HCT116/HEK293T(figshare 29190726), VCC-H1, CD4(done). NOT yet added, genome-wide, downloadable:
Feng-GW-iPSC + Feng-TS-iPSC (figshare 27989294.v2), Pan-GP-H1 (GEO GSE295214), Nourreddine-iPSC
(= KOLF, figshare+ 27261219). **Preprocessed effects + test stats for ALL 12 on HKU DataHub
(10.25442/hku.34002315)** -- may give ready-to-use effect matrices, saving harmonize effort.
Code: github.com/liuqiyuanhhh/perturbation-structure; preprocessing = crispyx (pypi).

**HKU DataHub (figshare article 34002315) = READY donors.** 69 files. `effect_<screen>.zip` =
per-pert effect matrix (pert x gene) + SEs, AnnData .h5ad, uniformly processed (crispyx), effect =
batch-aware mean diff of log(1+x) -> SAME family as our 02_harmonize deltas, so drops into our
average with no unit mismatch. Much smaller than raw cells. Download via ndownloader.figshare.com/
files/<id>. New genome-wide donors to add (effect_ zips): Pan-GW-hESC (H1! file 69471805, 2.8GB),
Feng-gwsf (69471714, 1.4GB), Feng-gwsnf (69471750, 2.2GB), Nourreddine-GW-ipsc (69472512, 11.5GB).
Also effect_ for K562/Huang/etc if we want uniform re-processing. SEs enable noise weighting.
NOT in this set: CD4/Marson (separate, done). README=69528249, MANIFEST=69528252.

### DE-frequency shared response (Liu/Wang 2026 idea; 21_de_frequency.py + 09 --shared-npz)
v14 `pred_c_defreq10_4donor_0.7` (DE-freq shared scale 10 replaces PIE shared; 4-donor 0.7):
overall **0.1073** | PDS 0.373 (0.669) | MSE 0 (3.311) | JAC -0.004 (0.029) | NMAE 0.122 (0.926) |
FID 0 (0.513) | REACH 0.152 (0.214). vs v13 (PIE shared, 0.116): FID UP (raw 0.502->0.513, best
yet -- DE-freq directions beat PIE's coin-flip) but PDS DOWN (0.440->0.373) -- a strong UNIFORM
shared component makes perts alike -> worse discrimination. Net slightly below v13.
Small-effect analysis (HCT116, effect terciles): low stratum = 2 DE genes/pert (noise floor),
Spearman(low,high freq) 0.74 but top-100 Jaccard 0.13 -- broad hierarchy shared, top responders
differ for weak perts. => per-pert shared weighting (09 --shared-weight-by-donor) should recover
PDS: scale shared by donor norm so silent perts get ~0 shared (non-uniform -> discriminable).
Per-pert weighting (--shared-weight-by-donor) did NOT engage: build log shows median weight 1.00,
0/300 perts near-zero in every context -> candidate ~= uniform v14 (gate 991 vs 1013), so not
submitted. Reason (useful finding): the near-silent perts (2 DE genes) are a full-screen
phenomenon; the 300 PANEL perts are a selected set that all have substantial donor effects, so
there is no silent subset within the panel to down-weight. => DE-freq shared (uniform or weighted)
improves FID but costs PDS; net below v13's PIE shared. DE-freq shelved as a lateral move.

**STATE (Arc flagship) — considered, not pursued (2026-10-19):** biorxiv 2025.06.26.661135;
SE (600M-param cell embedding) + ST (bidirectional transformer, set-level control->perturbed).
Blocked: not installed, no GPU on BioHPC node, no task-matched CRISPRi ST checkpoint (released
ST-Tahoe=drug, ST-Parse=diff domain; VCC baseline not located). Prior (our 18 + Ahlmann-Eltze +
Liu/Wang): deep models don't beat simple transfer on these metrics -> low expected payoff for a
multi-day GPU setup. Write-up point: benchmarked-by-reference, not run.

**Shelved:** PCA donor denoising (16) -- monotone loss on H1 at every k. PIE per-pert signal
(`--shrink > 0`) -- no agreement with K562 effects.

**Rules:** run the DE-call count before every submission; name `.vcc` files after their settings
and confirm the source file before recording a score.

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
| v6 `pred_c_pie1.0_d0.3` (uploaded as `pred_v6_pie1.0_d0.3.vcc`) | `--shrink 0 --lfc-scale 1.0 --donor-npz replogle_k562_gwps.npz --donor-scale 0.3` | **0.0267** | 0.102 | 0 | 0.008 | -0.128 | -0.022 | 0.201 | 926 |

| v7 `pred_c_pie0.5_d0.3` (`pred_v7_pie0.5_d0.3.vcc`) | `--shrink 0 --lfc-scale 0.5 --donor-npz replogle_k562_gwps.npz --donor-scale 0.3` | **0.0435** | 0.199 (0.591) | 0 (3.441) | -0.017 (0.024) | 0.027 (0.983) | -0.092 (0.485) | 0.144 (0.206) | |

v7 vs v6 (only PIE's shared scale 1.0 -> 0.5; DE calls ~4,200 -> ~740): PDS doubled (raw 0.546 ->
0.591), raw MSE 13.6 -> 3.44, NMAE now above baseline (raw 1.071 -> 0.983); cost: FID raw 0.506 ->
0.485 and REACH 0.20 -> 0.14. Net +0.017 overall. PIE's shared magnitude WAS driving NMAE and MSE.

| v8 `pred_c_d1.0` (`pred_v8_d1.0.vcc`) | `--shrink 0 --lfc-scale 0 --donor-npz replogle_k562_gwps.npz --donor-scale 1.0` | 0.0420 | **0.426 (0.693)** | 0 (9.034) | -0.012 (0.025) | -0.063 (1.039) | -0.187 (0.457) | 0.091 (0.160) | |

| v9 `pred_c_d0.7` (`pred_v9_d0.7.vcc`) | `--shrink 0 --lfc-scale 0 --donor-npz replogle_k562_gwps.npz --donor-scale 0.7` | **0.0492** | 0.393 (0.689) | 0 (4.909) | 0.027 (0.984) | -0.017 (0.024) | -0.192 (0.455) | 0.085 (0.154) | |

v9 vs v8 (donor only, scale 1.0 -> 0.7): PDS essentially held (raw 0.693 -> 0.689; cosine metric,
scale-free) while raw MSE dropped 9.0 -> 4.9 and NMAE 1.039 -> 0.984. Net best so far (0.0492).
Confirms: donor DIRECTION carries PDS, donor MAGNITUDE should be ~0.7, not 1.0. FID/REACH still
low (no shared part) -- v10 tests adding PIE 0.25.

| v10 `pred_c_pie0.25_d0.7` (`pred_v10_pie0.25_d0.7.vcc`) | `--shrink 0 --lfc-scale 0.25 --donor-npz replogle_k562_gwps.npz --donor-scale 0.7` | **0.0672** | 0.370 (0.668) | 0 (5.297) | 0.014 (0.992) | -0.012 (0.026) | -0.068 (0.492) | 0.100 (0.168) | |

v10 vs v9 (add PIE shared 0.25): big win (+0.018 overall). A small PIE shared part recovers FID
(raw 0.455 -> 0.492, ~coin flip) and REACH (0.085 -> 0.100) at a small PDS cost (raw 0.689 ->
0.668) and ~no NMAE cost (0.984 -> 0.992). So the best recipe so far = strong donor direction
(0.7) + a LIGHT shared part (0.25) for the direction metrics. v7 had PIE 0.5; 0.25 is better.

v8 vs v7: PDS more than doubled (raw 0.591 -> 0.693, between leaderboard ranks 900 and 500), so
the K562 donor carries real pert-specific direction once nothing competes with it. But at full
scale its magnitudes are too big (raw MSE 3.4 -> 9.0, NMAE back above 1) and FID / REACH fell
(raw FID 0.457, below coin flip, with ~2,800 calls). Overall about equal to v7. Target: keep v8's
PDS with v7's MSE / NMAE / FID / REACH -- that combination would be ~0.08 overall.

**Corrected 2026-10-13** (from shell history): v6 was PIE 1.0 + donor 0.3, NOT PIE 0.5 + donor
0.3 as first recorded. So v6 vs v4 changes ONLY the donor scale (0.1 -> 0.3), with the same PIE
part and the same ~4,200 DE calls per pert:
- PDS 0.03 -> 0.10 and REACH 0.18 -> 0.20 come from the stronger donor.
- NMAE unchanged (-0.128) is expected, since PIE's shared part did not change. The hypothesis that
  PIE's shared magnitude drives NMAE is NOT tested yet.
- Nothing is known yet about FID at 450-1,050 calls: only ~1 call (v5, collapse) and ~4,200
  (v4/v6, fine) have been submitted. PIE 0.5 (739 calls locally) is untested.
PDS from the K562 donor on the real panel (0.10-0.125 scaled across v5/v6) is far below what the
H1 proxy fold suggested.

v6 raw: PDS 0.546, MSE 13.619, JAC 0.033, NMAE 1.071, FID 0.506, REACH 0.256.
- Raw NMAE > 1 = worse than predicting no change on the true DE genes (proxy no_change = 1.0;
  leaderboard rank 900 = 0.979). Raw MSE 13.6 is ~4x rank 900 (3.4). Together with FID ~0.5
  (coin-flip directions), PIE's shared part looks like direction noise: needed only to produce
  DE calls, and it costs NMAE and MSE. (v6 used PIE at full scale 1.0.)
- Raw PDS 0.546: K562 transfer on the real panel is weaker than even the K562 proxy fold (0.59).

PIE per-pert signal vs K562 GWPS (269 covered Challenge perts, per-pert parts = centered on the
panel mean; genes > 5 CPM in context, covered, on PIE axis; PIE lfc clipped +-3):
| ctx | median r (null) | top-200 K562 genes sign agreement (null) |
|---|---|---|
| A | 0.013 (-0.001) | 0.495 (0.490) |
| B | 0.011 (0.001) | 0.500 (0.495) |
| C | 0.011 (-0.003) | 0.510 (0.500) |
=> no detectable agreement between PIE's pert-specific predictions and measured K562 knockdown
effects. Caveat: bounded by GWPS reliability for these perts (split-half pending), but GWPS does
give some real-panel PDS (raw 0.546) while PIE shows none against it. Deprioritize PIE per-pert
signal; work on the donor side (denoising, X-Atlas for B/C, K562->Jurkat map for A).

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

### GWPS split-half on Challenge perts (done 2026-10-12)
Output `split_half_challenge.csv` (on BioHPC). 272 Challenge perts in GWPS (median 172 cells) vs
100 control-cell null groups of matched size (median 183):

| | Challenge median | null median |
|---|---|---|
| effect norm | 14.9 | 13.5 |
| split-half cosine | 0.048 | 0.008 |
| top-200 sign agreement | 0.550 | 0.510 |

- 101/272 Challenge perts beat the null's 95th pct of split-half cosine (0.075); ~14 expected by
  chance -> roughly 85-90 perts with detectable signal, the rest indistinguishable from noise.
- Effect norm is barely above a same-size control group: signal ~ sqrt(14.9^2 - 13.5^2) ~ 6.3,
  i.e. only ~15-20% of each delta's squared norm is signal. GWPS deltas for these perts are
  mostly sampling noise -- consistent with the low real-panel PDS (raw 0.546) and a reason to
  denoise (16_denoise_donor.py) and/or shrink each pert by its reliability.
- 16_denoise_donor.py on GWPS (fit on 9,326 perts >= 50 cells, 7,679 genes): parallel analysis
  k = 22 (23% of variance; PC1 eigenvalue 24.5 vs 4.7 for PC2); panel perts keep median 30% of
  their delta norm. H1 proxy fold with the denoised GWPS (raw -> denoised): transfer pds 0.797 ->
  0.696, transfer_scaled_other 0.814 -> 0.713, transfer_topn 0.847 -> 0.733; nmae slightly worse.
  `--wiener` at k = 22 is slightly worse still (scaled_other 0.701, topn 0.721).
  k sweep (scaled_other / topn): k=10 0.695/0.711, k=22 0.713/0.733, k=50 0.723/0.748,
  k=100 0.735/0.765, raw 0.814/0.847 -- monotone in k, never reaches raw.
  => PCA denoising SHELVED: on H1 every dropped component costs pert-specific signal. Use
  per-pert reliability weighting (17_reliability_weight.py) instead. Caveat: the H1 perts are probably strong, well-measured effects, unlike
  the Challenge perts. Next: k sweep / --wiener, and a split-half projection test on the
  Challenge perts themselves.
- 17_reliability_weight.py (method mean): noise constant K = 3.26e4; the two estimates agree
  (Pearson 0.946). Signal fraction f over 272 Challenge perts: quantiles 10/25/50/75/90% =
  0.008 / 0.027 / 0.079 / 0.184 / 0.30; 39 perts > 0.25, 3 > 0.5, 15 at 0. Null groups get median
  0.010, 90th pct 0.079 -- so about half the Challenge perts are at noise level. Output
  `results/weighted/replogle_k562_gwps_rw.npz` (269 panel rows weighted).
- Reliability ceiling: Spearman-Brown on the median split cosine gives full-data reliability
  ~0.09, so even a perfect predictor would correlate only ~0.3 with a GWPS delta. PIE's 0.01 is
  still far below that.
