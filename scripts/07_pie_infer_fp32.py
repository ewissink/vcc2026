#!/usr/bin/env python
"""Run `pie infer` with the checkpoint's precision overridden to full float32 on CPU.

PIE reads trainer.precision from the checkpoint (bf16-mixed) and wraps the forward pass in
bfloat16 autocast. On CPUs without native bf16 that can be slow, so this wrapper patches the
loaded config before prediction. All other behaviour is PIE's own.

Usage (same arguments as `pie infer`):
    PIE_THREADS=16 python scripts/07_pie_infer_fp32.py experiment_name=replogle_xdataset \
        device=cpu rows_kind=query rows_path=$PWD/query.json \
        "preprocessed_dirs=[$PIE_DATA_ROOT/vcc2026_val]" output_path=$PWD/pred.parquet overwrite=true

Env: PIE_THREADS (default: all cores) sets torch/BLAS threads; PIE_PRECISION (default 32-true).
Not tested against the real pie package here; if the patch fails the error will say where.
"""
import os
import sys

n = os.environ.get("PIE_THREADS")
if n:
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[var] = n

import pie.predict as pp  # noqa: E402  (after thread env vars so BLAS picks them up)
import torch  # noqa: E402

if n:
    torch.set_num_threads(int(n))

PRECISION = os.environ.get("PIE_PRECISION", "32-true")
_orig_load = pp.load_checkpoint


def _load_checkpoint(path, *a, **kw):
    loaded = _orig_load(path, *a, **kw)
    trainer = loaded.config.trainer
    try:
        trainer.precision = PRECISION
    except Exception:  # frozen config object
        object.__setattr__(trainer, "precision", PRECISION)
    print(f"[07] precision overridden: {loaded.config.trainer.precision}; "
          f"torch threads: {torch.get_num_threads()}", flush=True)
    return loaded


pp.load_checkpoint = _load_checkpoint

from pie.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.argv = ["pie", "infer"] + sys.argv[1:]
    sys.exit(main())
