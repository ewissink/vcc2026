# CPU-only image for the VCC 2026 pipeline.
# Build:  docker build -t vcc2026 .
# Run:    docker run --rm -it --cpus=40 -v $(pwd)/data:/data vcc2026

FROM condaforge/miniforge3:latest

WORKDIR /workspace

COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

# Activate the env for all subsequent RUN/CMD instructions.
SHELL ["conda", "run", "-n", "vcc2026", "/bin/bash", "-c"]

# Thread env vars matching a 40-core host. Override at `docker run`
# time with -e if your host has a different core count.
ENV OMP_NUM_THREADS=40 \
    MKL_NUM_THREADS=40 \
    OPENBLAS_NUM_THREADS=40 \
    NUMEXPR_NUM_THREADS=40

COPY vcc2026_pipeline.py .

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "vcc2026", "python", "vcc2026_pipeline.py"]
