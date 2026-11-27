#!/usr/bin/env bash
# Restart ingest workers whose sync_cursors last_run_at is older than STALE_MIN.
# Respects PHI_OMICS_SNIFFER_STATE_DIR / OMICS_SNIFFER_DATA_DIR / OMICS_SNIFFER_LOGS_DIR.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
STATE_DIR="${PHI_OMICS_SNIFFER_STATE_DIR:-${OMICS_SNIFFER_STATE_DIR:-${OMICS_SNIFFER_ROOT:-$CODE_DIR}}}"
DATA_DIR="${OMICS_SNIFFER_DATA_DIR:-$STATE_DIR/data}"
LOGS_DIR="${OMICS_SNIFFER_LOGS_DIR:-$STATE_DIR/logs}"
export PHI_OMICS_SNIFFER_STATE_DIR="$STATE_DIR"

cd "$CODE_DIR"
mkdir -p "$DATA_DIR" "$LOGS_DIR"
STALE_MIN="${STALE_MIN:-15}"
LOG="$LOGS_DIR/ingest_hang_watchdog.log"
exec >>"$LOG" 2>&1
echo "$(date '+%F %T') hang-watchdog start stale=${STALE_MIN}m state=$STATE_DIR"

age_min() {
  local iso="$1"
  python3 - "$iso" <<'PY'
import sys
from datetime import datetime, timezone, timedelta
TZ = timezone(timedelta(hours=8))
s = sys.argv[1]
for fmt in ("%Y-%m-%dT%H:%M:%S%z",):
    try:
        t = datetime.strptime(s.replace("+08:00","+0800").replace("-08:00","-0800"), fmt)
        break
    except ValueError:
        t = None
else:
    print(99999); raise SystemExit
now = datetime.now(TZ)
print(int((now - t.astimezone(TZ)).total_seconds() // 60))
PY
}

cursor_finished() {
  local db="$1" source="$2"
  sqlite3 "$db" "SELECT cursor_value FROM sync_cursors WHERE source='$source';" 2>/dev/null | python3 -c '
import sys, json
raw = sys.stdin.read().strip()
if not raw:
    raise SystemExit(1)
try:
    c = json.loads(raw)
except Exception:
    raise SystemExit(1)
if c.get("done") is True or c.get("completed_at"):
    raise SystemExit(0)
raise SystemExit(1)
' >/dev/null 2>&1
}

check_one() {
  local name="$1" db="$2" source="$3" pattern="$4" start_cmd="$5"
  local iso age pid
  iso=$(sqlite3 "$db" "SELECT last_run_at FROM sync_cursors WHERE source='$source';" 2>/dev/null || true)
  if [[ -z "$iso" ]]; then
    echo "$(date '+%F %T') $name: no cursor yet; skip"
    return
  fi
  if cursor_finished "$db" "$source"; then
    echo "$(date '+%F %T') $name: cursor finished; skip restart"
    return
  fi
  age=$(age_min "$iso")
  pid=$(pgrep -f "$pattern" | head -1 || true)
  echo "$(date '+%F %T') $name: last_run_at=$iso age=${age}m pid=${pid:-none}"
  if (( age < STALE_MIN )); then
    return
  fi
  if [[ -n "$pid" ]]; then
    echo "$(date '+%F %T') $name: STALE — SIGTERM $pid"
    kill -TERM "$pid" 2>/dev/null || true
    for i in 1 2 3 4 5 6 7 8 9 10; do
      kill -0 "$pid" 2>/dev/null || break
      sleep 1
    done
    if kill -0 "$pid" 2>/dev/null; then
      echo "$(date '+%F %T') $name: SIGKILL $pid"
      kill -KILL "$pid" 2>/dev/null || true
      sleep 1
    fi
  fi
  if [[ "$name" == "AE" ]]; then
    if pgrep -f 'arrayexpress_full_scan.py' >/dev/null || pgrep -f 'arrayexpress_watchdog.sh' >/dev/null; then
      echo "$(date '+%F %T') $name: AE process/watchdog present; skip start"
      return
    fi
  fi
  if pgrep -f "$pattern" >/dev/null; then
    echo "$(date '+%F %T') $name: already running after kill; skip"
    return
  fi
  echo "$(date '+%F %T') $name: starting: $start_cmd"
  eval "$start_cmd"
}

while true; do
  check_one SRA "$DATA_DIR/omics_v1.db" sra_orphan_scan 'scanners/sra_orphan_scan.py' \
    "nohup python3 -u '$CODE_DIR/scanners/sra_orphan_scan.py' >> '$LOGS_DIR/sra_orphan_scan.log' 2>&1 &"
  check_one AE "$DATA_DIR/omics_ae_staging.db" arrayexpress_full_scan 'scanners/arrayexpress_full_scan.py' \
    "nohup bash '$CODE_DIR/scripts/arrayexpress_watchdog.sh' &"
  check_one PRIDE "$DATA_DIR/omics_pride_staging.db" pride_full_scan 'scanners/pride_full_scan.py' \
    "nohup python3 -u '$CODE_DIR/scanners/pride_full_scan.py' >> '$LOGS_DIR/pride_full_scan.log' 2>&1 &"
  check_one MTBLS "$DATA_DIR/omics_metabolights_staging.db" metabolights_full_scan 'scanners/metabolights_full_scan.py' \
    "nohup python3 -u '$CODE_DIR/scanners/metabolights_full_scan.py' >> '$LOGS_DIR/metabolights_full_scan.log' 2>&1 &"
  check_one MW "$DATA_DIR/omics_mw_staging.db" mw_full_scan 'scanners/mw_full_scan.py' \
    "nohup python3 -u '$CODE_DIR/scanners/mw_full_scan.py' >> '$LOGS_DIR/mw_full_scan.log' 2>&1 &"
  check_one MASSIVE "$DATA_DIR/omics_massive_staging.db" massive_full_scan 'scanners/massive_full_scan.py' \
    "nohup python3 -u '$CODE_DIR/scanners/massive_full_scan.py' >> '$LOGS_DIR/massive_full_scan.log' 2>&1 &"
  check_one CXG "$DATA_DIR/omics_cellxgene_staging.db" cellxgene_full_scan 'scanners/cellxgene_full_scan.py' \
    "nohup python3 -u '$CODE_DIR/scanners/cellxgene_full_scan.py' >> '$LOGS_DIR/cellxgene_full_scan.log' 2>&1 &"
  check_one GDC "$DATA_DIR/omics_gdc_tcga_staging.db" gdc_tcga_sample_scan 'scanners/gdc_tcga_sample_scan.py' \
    "nohup python3 -u '$CODE_DIR/scanners/gdc_tcga_sample_scan.py' >> '$LOGS_DIR/gdc_tcga_sample_scan.log' 2>&1 &"
  check_one GTEX "$DATA_DIR/omics_gtex_staging.db" gtex_sample_scan 'scanners/gtex_sample_scan.py' \
    "nohup python3 -u '$CODE_DIR/scanners/gtex_sample_scan.py' >> '$LOGS_DIR/gtex_sample_scan.log' 2>&1 &"
  check_one CPTAC "$DATA_DIR/omics_cptac_staging.db" cptac_sample_scan 'scanners/cptac_sample_scan.py' \
    "nohup python3 -u '$CODE_DIR/scanners/cptac_sample_scan.py' >> '$LOGS_DIR/cptac_sample_scan.log' 2>&1 &"
  check_one ICPC "$DATA_DIR/omics_icpc_staging.db" icpc_sample_scan 'scanners/icpc_sample_scan.py' \
    "nohup python3 -u '$CODE_DIR/scanners/icpc_sample_scan.py' >> '$LOGS_DIR/icpc_sample_scan.log' 2>&1 &"
  sleep 300
done
