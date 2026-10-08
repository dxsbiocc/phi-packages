# syntax=docker/dockerfile:1
# Shared R image for every wrapper under
# resources/wrappers/modules/local/differential-expression/. Built and
# tagged locally (not pushed to any registry) -- see environment.yml
# beside this file for exactly what's in it and why.
#
# Build (from this directory):
#   docker build -t phi/differential-expression-r:1.0.0 .
#
# Point a wrapper's ../main.nf `container` directive at this same tag to
# switch it over.
FROM condaforge/miniforge3:24.11.3-2

COPY environment.yml /opt/environment.yml
# A separate named env, not `-n base` -- updating base in place upgrades
# libraries mamba's own CLI is linked against (broke libxml2 the first
# time, crashing `mamba clean` right after a successful install).
# The cache mount on /opt/conda/pkgs persists already-downloaded packages
# across retries of this RUN -- this network has dropped mid-download
# (IncompleteRead) on every attempt so far, and without this a retry
# re-downloads everything from scratch instead of resuming.
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -p /opt/conda/envs/de -f /opt/environment.yml \
    && rm /opt/environment.yml

ENV PATH=/opt/conda/envs/de/bin:$PATH
