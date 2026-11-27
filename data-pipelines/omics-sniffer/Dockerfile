FROM python:3.12-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        sqlite3 \
        util-linux \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/phi-omics-sniffer
COPY pyproject.toml schema_v1.sql ./
COPY pipeline ./pipeline
COPY scanners ./scanners
COPY scripts ./scripts

RUN pip install --no-cache-dir -e . \
    && chmod +x scripts/*.sh \
    && mkdir -p /var/lib/phi-omics-sniffer/serving /var/lib/phi-omics-sniffer/data /var/lib/phi-omics-sniffer/logs

ENV PHI_OMICS_SNIFFER_STATE_DIR=/var/lib/phi-omics-sniffer \
    PHI_OMICS_SNIFFER_DATA_DIR=/var/lib/phi-omics-sniffer/serving \
    PYTHONUNBUFFERED=1

ENTRYPOINT ["/bin/bash", "/opt/phi-omics-sniffer/scripts/run_sync_pipeline.sh"]
