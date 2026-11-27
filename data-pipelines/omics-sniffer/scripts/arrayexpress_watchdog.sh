#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
STATE_DIR="${PHI_OMICS_SNIFFER_STATE_DIR:-${OMICS_SNIFFER_STATE_DIR:-${OMICS_SNIFFER_ROOT:-$CODE_DIR}}}"
LOGS_DIR="${OMICS_SNIFFER_LOGS_DIR:-$STATE_DIR/logs}"
export PHI_OMICS_SNIFFER_STATE_DIR="$STATE_DIR"

cd "$CODE_DIR"
mkdir -p "$LOGS_DIR"
exec >> "$LOGS_DIR/arrayexpress_full_scan.log" 2>&1
echo "$(date '+%F %T') watchdog start (state=$STATE_DIR)"
while true; do
  python3 -u "$CODE_DIR/scanners/arrayexpress_full_scan.py"
  code=$?
  echo "$(date '+%F %T') scan exited code=$code"
  if [[ $code -eq 0 ]]; then
    echo "$(date '+%F %T') scan completed; watchdog stop"
    exit 0
  fi
  echo "$(date '+%F %T') retry in 120s"
  sleep 120
done
