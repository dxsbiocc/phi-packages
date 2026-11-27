"""Shared helpers for metadata-only staging scanners (iProX/jPOST, DepMap, ICGC).
Never touches omics_v1.db."""
from __future__ import annotations
import hashlib, json, logging, re, sqlite3, sys, time
from datetime import datetime, timedelta, timezone
from pathlib import Path
try:
    from .http_util import TIMEOUT, http_get
    from .paths import CODE_ROOT, SCHEMA_PATH, data_dir, logs_dir
except ImportError:
    from http_util import TIMEOUT, http_get
    from paths import CODE_ROOT, SCHEMA_PATH, data_dir, logs_dir

ROOT = CODE_ROOT
DATA_DIR = data_dir()
LOGS_DIR = logs_dir()
SCHEMA = SCHEMA_PATH
TZ = timezone(timedelta(hours=8))
MAX_RETRIES = 8
LIVER_RE = re.compile(r"liver|hepat|\bhcc\b|cholangi|bile duct|biliary|hepatoblastoma|gallbladder", re.I)

def now_iso(): return datetime.now(TZ).strftime("%Y-%m-%dT%H:%M:%S%z")
def chash(o): return hashlib.sha256(json.dumps(o, sort_keys=True, default=str).encode()).hexdigest()[:16]
def tnorm(t): return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9 ]", " ", (t or "").lower())).strip()
def is_liver(*p): return bool(LIVER_RE.search(" ".join(str(x) for x in p if x)))

def get_logger(name):
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(LOGS_DIR / f"{name}.log")])
    return logging.getLogger(name)

class Throttle:
    def __init__(s, gap): s.gap, s.t = gap, 0.0
    def __call__(s):
        d = time.monotonic() - s.t
        if d < s.gap: time.sleep(s.gap - d)
        s.t = time.monotonic()

def fetch(url, log, thr, params=None, ok404=False, **kw):
    for a in range(1, MAX_RETRIES + 1):
        thr()
        try:
            r = http_get(url, params=params, timeout=TIMEOUT, **kw)
            if r.status_code == 404 and ok404: return None
            if r.status_code in (429,) or r.status_code >= 500: raise RuntimeError(f"HTTP {r.status_code}")
            r.raise_for_status(); return r
        except Exception as e:
            log.warning("retry %d/%d %s: %s", a, MAX_RETRIES, url, e); time.sleep(min(60, 2 ** a))
    raise RuntimeError(f"giving up {url}")

def connect(db: Path, extra_ddl: str):
    assert db.name != "omics_v1.db"
    db.parent.mkdir(parents=True, exist_ok=True)
    new = not db.exists()
    c = sqlite3.connect(db, timeout=30); c.execute("PRAGMA journal_mode=WAL"); c.execute("PRAGMA foreign_keys=ON")
    if new: c.executescript(SCHEMA.read_text())
    c.executescript(extra_ddl); c.commit(); return c

def load_cursor(c, src):
    r = c.execute("SELECT cursor_value FROM sync_cursors WHERE source=?", (src,)).fetchone()
    return json.loads(r[0]) if r and r[0] else {}

def save_cursor(c, src, cur, err=None):
    c.execute("INSERT INTO sync_cursors(source,cursor_value,cursor_kind,last_run_at,last_error) VALUES(?,?,?,?,?) "
              "ON CONFLICT(source) DO UPDATE SET cursor_value=excluded.cursor_value,last_run_at=excluded.last_run_at,last_error=excluded.last_error",
              (src, json.dumps(cur), "json", now_iso(), err)); c.commit()

def upsert_dataset(c, d: dict):
    ts = now_iso(); d = {**d}; d.setdefault("title_norm", tnorm(d.get("title")))
    d["content_hash"] = chash(d); d.update(updated_at=ts, last_checked=ts)
    for k in ("disease_tags", "platforms"):
        if isinstance(d.get(k), (list, dict)): d[k] = json.dumps(d[k], ensure_ascii=False)
    cols = list(d)
    c.execute(f"INSERT INTO datasets({','.join(cols)},first_seen_at) VALUES({','.join('?'*len(cols))},?) "
              f"ON CONFLICT(canonical_id) DO UPDATE SET {','.join(f'{k}=excluded.{k}' for k in cols)}",
              [d[k] for k in cols] + [ts])

def upsert_sample(c, s: dict):
    s = {**s}
    if isinstance(s.get("attributes"), dict): s["attributes"] = json.dumps(s["attributes"], ensure_ascii=False)
    s["content_hash"] = chash(s); s["updated_at"] = now_iso(); cols = list(s)
    c.execute(f"INSERT OR REPLACE INTO samples({','.join(cols)}) VALUES({','.join('?'*len(cols))})", [s[k] for k in cols])

def add_member(c, ds, source, acc, role, url=None):
    c.execute("INSERT OR IGNORE INTO dataset_members(dataset_id,source,accession,role,external_url) VALUES(?,?,?,?,?)", (ds, source, acc, role, url))

def add_artifact(c, a: dict):
    a = {**a, "last_checked": now_iso()}; cols = list(a)
    c.execute(f"INSERT OR REPLACE INTO dataset_artifacts({','.join(cols)}) VALUES({','.join('?'*len(cols))})", [a[k] for k in cols])
