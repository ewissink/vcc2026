#!/usr/bin/env python
"""Run PIE inference for a query file as several parallel processes and merge the result.

Splits {"<dataset>.<context>": [perts]} into N chunks of (context, perturbation) pairs, runs one
scripts/07_pie_infer_fp32.py process per chunk (each pays the ~2 min start-up and loads its own
copy of the model/evidence, so watch memory), then concatenates the parquet files with pyarrow so
the gene axis stored in the file metadata is kept. Chunks whose output exists are skipped, so
re-running after a failure only redoes the missing ones.

Usage:
    python scripts/08_pie_infer_parallel.py --query query.json --out pred.parquet \
        --workdir pie_chunks --n-chunks 6 --threads 8 \
        --pie-arg experiment_name=replogle_xdataset --pie-arg "preprocessed_dirs=[$PIE_DATA_ROOT/vcc2026_val]"
(--pie-arg is repeated; each value is passed to `pie infer` as given. device=cpu, rows_kind=query,
rows_path, output_path and overwrite are set here.)
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def split_query(query, n):
    pairs = [(k, p) for k, perts in query.items() for p in perts]
    chunks = [pairs[i::n] for i in range(n)]          # round-robin keeps contexts balanced
    out = []
    for c in chunks:
        q = {}
        for k, p in c:
            q.setdefault(k, []).append(p)
        out.append(q)
    return [q for q in out if q], len(pairs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--n-chunks", type=int, default=6)
    ap.add_argument("--threads", type=int, default=8, help="PIE_THREADS per process")
    ap.add_argument("--pie-arg", action="append", default=[])
    ap.add_argument("--worker", default=str(Path(__file__).with_name("07_pie_infer_fp32.py")))
    args = ap.parse_args()

    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    query = json.load(open(args.query))
    chunks, total = split_query(query, args.n_chunks)
    print(f"{total} (context, perturbation) pairs in {len(chunks)} chunks", flush=True)

    procs = []
    for i, q in enumerate(chunks):
        qpath, opath = work / f"chunk{i}.json", work / f"chunk{i}.parquet"
        json.dump(q, open(qpath, "w"))
        if opath.exists():
            print(f"chunk {i}: output exists, skipping", flush=True)
            continue
        cmd = [sys.executable, args.worker, *args.pie_arg, "device=cpu", "rows_kind=query",
               f"rows_path={qpath}", f"output_path={opath}", "overwrite=true"]
        log = open(work / f"chunk{i}.log", "w")
        env = dict(os.environ, PIE_THREADS=str(args.threads))
        procs.append((i, subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)))
        print(f"chunk {i}: started ({sum(len(v) for v in q.values())} pairs)", flush=True)

    failed = [i for i, p in procs if p.wait() != 0]
    if failed:
        raise SystemExit(f"chunks failed: {failed}; see {work}/chunk<i>.log, fix, and re-run")

    tables = [pq.read_table(work / f"chunk{i}.parquet") for i in range(len(chunks))]
    merged = pa.concat_tables(tables)
    keys = list(zip(merged["context"].to_pylist(), merged["perturbation"].to_pylist()))
    if merged.num_rows != total or len(set(keys)) != total:
        raise SystemExit(f"merge check failed: {merged.num_rows} rows, {len(set(keys))} unique, "
                         f"expected {total}")
    pq.write_table(merged, args.out)
    print(f"wrote {merged.num_rows} rows -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
