#!/usr/bin/env bash
# Final-entry builder: run the v17 recipe (4-donor average + PIE shared) across a donor/lfc scale
# grid on the final-test contexts, gate each by DE-call count, and print a summary to choose from.
# It does NOT submit (that's deliberate) and does NOT run PIE (different env) -- see CHECKLIST_OCT22.md.
#
# Prerequisites (produced by the checklist's upstream steps):
#   $DATA_DIR/gene_names.csv, $DATA_DIR/pert_counts.csv, $DATA_DIR/context_<NAME>.h5ad (new test controls)
#   $PRED  = PIE predictions for the new contexts (parquet)
#   $PIE_META = meta.json of the new controls-only PIE prep (gene axis)
#   the 4 donor .npz (unchanged from the leaderboard phase)
#
# Usage: edit the CONFIG block, then:  bash scripts/run_final.sh
set -euo pipefail

############################## CONFIG (edit for the final data) ##############################
DATA_DIR="data_dir_final"                               # context_*.h5ad, gene_names.csv, pert_counts.csv
PRED="pred_final_full.parquet"                          # PIE predictions on the new contexts
PIE_META="/workdir/emw97/pie/data/vcc2026_final/meta.json"
OUT="final_candidates"
WORKERS=32
DONORS="results/harmonized/replogle_k562_gwps.npz,results/harmonized/xatlas_hct116.npz,results/harmonized/xatlas_hek293t.npz,results/harmonized/cd4_rest.npz"
DONOR_SCALES=(0.9 1.0 1.1)                               # sweep (v17 best was 0.9; PDS still rising)
LFC_SCALES=(0.25 0.20)                                   # shared scale (0.25 = v17)
MIN_CALLS=500                                            # DE-call floor per context (FID yield guard)
##############################################################################################

mkdir -p "$OUT"
echo "== preflight =="
for f in "$DATA_DIR/gene_names.csv" "$DATA_DIR/pert_counts.csv" "$PRED" "$PIE_META"; do
  [ -e "$f" ] || { echo "MISSING: $f"; exit 1; }
done
ls "$DATA_DIR"/context_*.h5ad >/dev/null || { echo "MISSING: $DATA_DIR/context_*.h5ad"; exit 1; }
IFS=',' read -ra DARR <<< "$DONORS"
for d in "${DARR[@]}"; do [ -e "$d" ] || { echo "MISSING donor: $d"; exit 1; }; done
# gene axis must match across gene_names, donors, and the PIE meta
python - "$DATA_DIR/gene_names.csv" "$PIE_META" "${DARR[@]}" <<'PY'
import sys, json, numpy as np
gn = sys.argv[1]; meta = sys.argv[2]; donors = sys.argv[3:]
genes = [l.strip().split(",")[0] for l in open(gn) if l.strip()]
if genes[0].lower() in ("gene","gene_name","gene_names","symbol"): genes = genes[1:]
G = set(genes); print(f"gene_names: {len(genes)} genes")
m = json.load(open(meta)); pie = set(m["genes"]); print(f"PIE axis: {len(pie)} genes; {len(pie & G)} in gene_names")
for d in donors:
    z = np.load(d, allow_pickle=True)
    ok = list(z["genes"].astype(str)) == genes
    print(f"donor {d.split('/')[-1]}: gene axis {'MATCHES' if ok else 'DIFFERS!!'}")
    assert ok, f"{d} gene axis != gene_names -- rebuild donor on the final gene_names order"
print("preflight OK")
PY

echo "== build scale grid =="
BUILT=()
for d in "${DONOR_SCALES[@]}"; do for l in "${LFC_SCALES[@]}"; do
  out="$OUT/pred_d${d}_l${l}.h5ad"
  echo "--- donor $d  lfc $l -> $out"
  python scripts/09_make_cells.py --pred "$PRED" --pie-meta "$PIE_META" \
    --data-dir "$DATA_DIR" --pert-counts "$DATA_DIR/pert_counts.csv" --out "$out" \
    --workers "$WORKERS" --shrink 0 --lfc-scale "$l" --donor "*=${DONORS}:${d}"
  BUILT+=("$out")
done; done

echo "== DE-call gate (per-context medians; floor=$MIN_CALLS) =="
python - "$DATA_DIR" "$MIN_CALLS" "${BUILT[@]}" <<'PY'
import sys, os, anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import mannwhitneyu
data_dir, floor = sys.argv[1], float(sys.argv[2]); files = sys.argv[3:]
rng = np.random.default_rng(0)
def lognorm(X):
    X = sp.csr_matrix(X).astype(np.float64)
    X = sp.diags(1e4/np.maximum(np.asarray(X.sum(1)).ravel(),1))@X; X.data=np.log1p(X.data); return X
def bh(p):
    n=len(p); o=np.argsort(p); q=p[o]*n/np.arange(1,n+1)
    q=np.minimum.accumulate(q[::-1])[::-1]; out=np.empty(n); out[o]=np.minimum(q,1); return out
ctrl={}
for f in files:
    a=ad.read_h5ad(f,backed="r"); obs=a.obs[["context","target_gene"]].astype(str); res=[]
    for ctx in sorted(obs.context.unique()):
        if ctx not in ctrl:
            c=ad.read_h5ad(f"{data_dir}/context_{ctx}.h5ad")
            Xc=lognorm(c.X[rng.choice(c.n_obs,min(2000,c.n_obs),replace=False)])
            cpm=np.asarray(Xc.expm1().mean(0)).ravel()*100; g=np.where(cpm>5)[0]; ctrl[ctx]=(Xc[:,g].toarray(),g)
        C,g=ctrl[ctx]
        perts=np.random.default_rng(1).choice(obs.target_gene[obs.context==ctx].unique(),25,replace=False)
        for p in perts:
            pos=np.sort(np.where((obs.context==ctx).values&(obs.target_gene==p).values)[0])
            P=lognorm(a[pos].to_memory().X)[:,g].toarray()
            res.append((ctx,int((bh(mannwhitneyu(P,C,axis=0).pvalue)<0.05).sum())))
    df=pd.DataFrame(res,columns=["ctx","n_sig"]); bycx=df.groupby("ctx").n_sig.median().to_dict()
    ok="PASS" if min(bycx.values())>=floor else "LOW"
    print(f"{os.path.basename(f)}: median {df.n_sig.median():.0f}  by ctx {bycx}  [{ok}]", flush=True)
PY

cat <<EOF

== next steps (manual) ==
1. Pick the scale by the LOCAL HARNESS score (cell-eval2), not just calls:
   for each candidate, score with scripts/15_score_local.py on harness_K562gw (proxy, has truth).
   Prefer the highest harness score among candidates marked PASS above.
2. Submit the chosen one:
   vcc prep $OUT/pred_d<D>_l<L>.h5ad -g $DATA_DIR/gene_names.csv --perts $DATA_DIR/pert_counts.csv -o final.vcc
   vcc submit final.vcc -m final_d<D>_l<L>
   (Only your LAST final entry counts -- submit the chosen one last.)
EOF
