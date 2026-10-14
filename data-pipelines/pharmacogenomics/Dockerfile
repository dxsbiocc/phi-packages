# R dependency for PharmacoSet RDS files.
#
# Build from this directory:
#   docker build -t phi-pgx:batch2 .
# Run an extract without the Python pipeline:
#   docker run --rm -v "$PWD/raw:/data" phi-pgx:batch2 \
#     Rscript /opt/phi-pgx/r/extract_pset.R /data/orcestra_gcsi/2019/gCSI2.rds /data/orcestra_gcsi/2019/extract
#
# Python entrypoint:
#   docker run --rm -v "$PWD:/work" -w /work phi-pgx:batch2 process
FROM rocker/r-ver:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-venv \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && R -e "BiocManager::install('PharmacoGx', ask=FALSE, update=FALSE)"

WORKDIR /opt/phi-pgx
COPY pyproject.toml README.md sources.yaml ./
COPY src ./src
COPY r ./r
ENV PYTHONPATH=/opt/phi-pgx/src
RUN pip3 install --break-system-packages duckdb pyyaml

ENTRYPOINT ["python3", "-m", "pgx.cli"]
