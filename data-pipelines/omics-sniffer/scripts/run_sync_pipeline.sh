#!/usr/bin/env bash
# Idempotent Linux server sync pipeline (Parquet-native by default, code/data separated):
#   Default: Run tiered incremental Parquet sync (no 15 GB main SQLite required)
#   --legacy-sqlite: Run legacy Staging SQLite -> omics_v1.db -> Parquet export
#
# Environment variables:
#   PHI_OMICS_SNIFFER_STATE_DIR / OMICS_SNIFFER_STATE_DIR: mutable runtime root (data/, serving/, logs/)
#   PHI_OMICS_SNIFFER_DATA_DIR / OMICS_SNIFFER_SERVING_DIR: Parquet serving root shared with phi-omics-gateway
#   OMICS_SNIFFER_LOGS_DIR: logs & lockfile directory
#
# Uses flock (when available) or mkdir lock to guarantee at most one instance runs at a time.
set -euo pipefail

MODE="parquet"
EXTRA_ARGS=()
for arg in "$@"; do
  case "$arg" in
    --legacy-sqlite) MODE="sqlite" ;;
    *) EXTRA_ARGS+=("$arg") ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
STATE_DIR="${PHI_OMICS_SNIFFER_STATE_DIR:-${OMICS_SNIFFER_STATE_DIR:-${OMICS_SNIFFER_ROOT:-$CODE_DIR}}}"
LOGS_DIR="${OMICS_SNIFFER_LOGS_DIR:-$STATE_DIR/logs}"
export PHI_OMICS_SNIFFER_STATE_DIR="$STATE_DIR"

cd "$CODE_DIR"
mkdir -p "$LOGS_DIR"
LOCK_FILE="$LOGS_DIR/sync_pipeline.lock"
LOG_FILE="$LOGS_DIR/sync_pipeline.log"

if [[ -n "${PYTHON_BIN:-}" && -x "${PYTHON_BIN}" ]]; then
  PY="$PYTHON_BIN"
elif [[ -x "$STATE_DIR/.venv/bin/python" ]]; then
  PY="$STATE_DIR/.venv/bin/python"
elif [[ -x "$CODE_DIR/.venv/bin/python" ]]; then
  PY="$CODE_DIR/.venv/bin/python"
else
  PY="$(command -v python3)"
fi

if command -v flock >/dev/null 2>&1; then
  exec 9>"$LOCK_FILE"
  if ! flock -n 9; then
    echo "$(date '+%F %T') [skip] another sync_pipeline instance is running" >> "$LOG_FILE"
    exit 0
  fi
fi

{
  echo "=========================================="
  echo "$(date '+%F %T') [start] sync_pipeline (code=$CODE_DIR, state=$STATE_DIR, mode=$MODE)"

  if [[ "$MODE" == "parquet" ]]; then
    "$PY" -u "$CODE_DIR/pipeline/incremental_parquet_sync.py" "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"
  else
    "$PY" -u "$CODE_DIR/pipeline/merge_staging_to_main.py" --all-finished --batch-size 1000
    "$PY" -u "$CODE_DIR/pipeline/export_parquet.py"
  fi

  echo "$(date '+%F %T') [done] sync_pipeline finished successfully"
} >> "$LOG_FILE" 2>&1
