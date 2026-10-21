# Final-entry checklist (data drops 2026-10-22, deadline 2026-11-05, NO leaderboard in between)

Best method = **v17**: 4-donor average (K562 GWPS, X-Atlas HCT116, X-Atlas HEK293T, CD4-rest),
PIE shared at `--lfc-scale 0.25`, `--donor-scale 0.9`, `--shrink 0`. Leaderboard score 0.1297.
Donors are REUSED unchanged (public data). Only the contexts + perts + PIE predictions are new.

## 0. Before the data drops (do now)
- [ ] `git pull` on BioHPC so scripts 09/17/19/20/run_final are current.
- [ ] Confirm the 4 donor `.npz` exist in `results/harmonized/` (+ `cd4_rest.npz`).
- [ ] Confirm the K562 harness exists (`harness_K562gw/` real.h5ad etc.) for blind scale selection.
- [ ] Confirm the PIE env + checkpoint still work (`replogle_xdataset/best_auprc.ckpt`).

## 1. Get + inspect the final data
- [ ] Download the final-test control profiles (3 new cell lines) + the 300 final perts + gene list.
- [ ] Check gene list == our 18,533 `gene_names.csv` order. If different, donors must be rebuilt on
      the new order (02/19/20) -- do that first.
- [ ] Check whether the 300 final perts == the val perts. If NEW, note donor panel coverage will
      differ (re-check with the coverage one-liner); if SAME, coverage already known.
- [ ] Build `data_dir_final/`: `context_<NAME>.h5ad` (the new controls, one per test line),
      `gene_names.csv`, `pert_counts.csv` (the 300 final perts). Match obs/var format of the val data.

## 2. Re-run PIE inference on the new contexts  (pie env -- the slow step, budget time)
- [ ] `pie prep` controls-only on the new control cells -> new prepped dir + meta.json
      (same command you used for vcc2026_val; set PIE_DATA_ROOT to a new vcc2026_final dir).
- [ ] `python scripts/06_pie_gene_filter.py` if the gene set needs the PIE-text filter.
- [ ] `python scripts/05_make_pie_query.py --pert-counts data_dir_final/pert_counts.csv --meta <new meta.json> --out query_final.json`
- [ ] `python scripts/08_pie_infer_parallel.py --query query_final.json --out pred_final_full.parquet ...`
      (CPU + parallel; this is the bottleneck.)

## 3. Build candidates + gate  (main env)
- [ ] Edit the CONFIG block of `scripts/run_final.sh` (DATA_DIR, PRED, PIE_META paths).
- [ ] `bash scripts/run_final.sh` -> builds the donor{0.9,1.0,1.1} x lfc{0.25,0.20} grid into
      `final_candidates/` and prints per-context DE-call medians with PASS/LOW.
- [ ] Discard any candidate marked LOW (a context below 500 calls -> FID yield risk).

## 4. Pick the scale BLIND via the local harness  (cell-eval2 env)
- [ ] For each PASS candidate's *scale*, build the matching cells on the K562 harness context and
      score with `scripts/15_score_local.py` against `harness_K562gw/real.h5ad` (append controls
      with 14 first). The harness is a different cell line, so use it to rank SCALES, not for an
      absolute number.
- [ ] Choose the scale with the best harness score among PASS candidates. (Prior: donor 0.9-1.0;
      PDS was still rising at 0.9 on the leaderboard.)

## 5. Submit  (only the LAST final entry counts)
- [ ] `vcc prep final_candidates/pred_d<D>_l<L>.h5ad -g data_dir_final/gene_names.csv --perts data_dir_final/pert_counts.csv -o final.vcc`
- [ ] `vcc submit final.vcc -m final_d<D>_l<L>`
- [ ] Submit the CHOSEN candidate LAST (earlier final entries are overwritten for scoring).
- [ ] Record the final choice + settings in EXPERIMENTS.md.

## Gotchas
- No leaderboard feedback Oct22-Nov5 -> the harness + DE-call gate are your only signal; don't guess.
- PIE inference time is the schedule risk -- start step 2 the moment data lands.
- 2 final entries/day max; only the last counts. Don't burn the slot before the harness pick.
- If a context's controls are tiny or oddly formatted, 09's expressed-gene filter and resampling
  still work, but sanity-check `context_*.h5ad` loads and has raw counts.
