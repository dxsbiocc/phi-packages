"""Shared contract logic for the deseq2-qc skill (stdlib only)."""

import hashlib
import io
import json
import math
import os
import re
import sys
import tempfile

QC_COLUMNS = [
    "schemaVersion",
    "disposition",
    "contrast",
    "geneIdColumn",
    "padjThreshold",
    "log2fcThreshold",
    "evaluated",
    "hits",
    "resultMd5",
]
DISPOSITIONS = ("eligible", "skip", "invalid")
MAX_LABEL = 200
MAX_REASONS = 10
UNNAMED_FIRST = "(unnamed-column-1)"
_NUM_RE = re.compile(r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$")
_MD5_RE = re.compile(r"^[0-9a-f]{32}$")


def emit(obj, code=0):
    sys.stdout.write(json.dumps(obj, separators=(",", ":"), ensure_ascii=False) + "\n")
    sys.exit(code)


def needs_input(reasons, question=None):
    obj = {"status": "needs_input", "reasons": reasons}
    if question is not None:
        obj["question"] = question
    emit(obj, 0)


def display_path(path):
    try:
        rel = os.path.relpath(path)
    except ValueError:
        return path
    return path if rel.startswith("..") else rel


def parse_cell(text):
    """Return float, or None for NA. Raise ValueError for invalid cells."""
    if text == "NA":
        return None
    if not _NUM_RE.match(text):
        raise ValueError(text)
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(text)
    return value


def _invalid(reasons, label="(unknown)"):
    return {
        "disposition": "invalid",
        "evaluated": 0,
        "hits": 0,
        "geneIdColumn": label,
        "hit_ids": set(),
        "evaluable_ids": set(),
        "reasons": reasons[:MAX_REASONS],
    }


def evaluate_results(data, gene_id_col, padj_threshold, log2fc_threshold):
    """Recompute evaluated/hit counts against the exact results bytes."""
    return evaluate_lines(io.BytesIO(data), gene_id_col, padj_threshold, log2fc_threshold)


def evaluate_lines(lines, gene_id_col, padj_threshold, log2fc_threshold):
    try:
        header = next(lines, b"").decode("utf-8").rstrip("\r\n").split("\t")
    except UnicodeDecodeError:
        return _invalid(["results file is not valid UTF-8"])
    if len(header) < 2:
        return _invalid(["header must have at least 2 tab-separated columns"])
    if any(name == "" for name in header[1:]):
        return _invalid(["only the first header column may be unnamed"])
    if len(set(header)) != len(header):
        return _invalid(["duplicate column names in header"])
    if gene_id_col:
        if gene_id_col not in header:
            return _invalid(["declared gene id column not found in header"])
        gi = header.index(gene_id_col)
        label = gene_id_col
    else:
        gi = 0
        label = header[0] or UNNAMED_FIRST
    for required in ("log2FoldChange", "padj"):
        if required not in header:
            return _invalid(["required column %r missing from header" % required], label)
    li = header.index("log2FoldChange")
    pi = header.index("padj")
    reasons = []
    seen = set()
    evaluable = set()
    hit_ids = set()
    evaluated = 0
    hits = 0
    for lineno, raw in enumerate(lines, start=2):
        if len(reasons) >= MAX_REASONS:
            return _invalid(reasons, label)
        try:
            line = raw.decode("utf-8").rstrip("\r\n")
        except UnicodeDecodeError:
            return _invalid(["results file is not valid UTF-8"], label)
        fields = line.split("\t")
        if len(fields) != len(header):
            reasons.append("row %d: expected %d fields, got %d" % (lineno, len(header), len(fields)))
            continue
        gid = fields[gi].strip()
        if gid == "" or gid == "NA":
            reasons.append("row %d: missing gene id" % lineno)
            continue
        if gid in seen:
            reasons.append("row %d: duplicate gene id" % lineno)
            continue
        seen.add(gid)
        try:
            lfc = parse_cell(fields[li])
            padj = parse_cell(fields[pi])
        except ValueError:
            reasons.append("row %d: invalid numeric cell (NA is the only valid missing)" % lineno)
            continue
        if padj is not None and not (0.0 <= padj <= 1.0):
            reasons.append("row %d: padj outside [0,1]" % lineno)
            continue
        if lfc is None or padj is None:
            continue
        evaluated += 1
        evaluable.add(gid)
        if padj < padj_threshold and abs(lfc) >= log2fc_threshold:
            hits += 1
            hit_ids.add(gid)
    if reasons:
        return _invalid(reasons, label)
    if not seen:
        return _invalid(["results contain no data records"], label)
    if evaluated == 0:
        return _invalid(["no evaluable genes: all padj/log2FoldChange values are NA"], label)
    disposition = "eligible" if hits > 0 else "skip"
    return {
        "disposition": disposition,
        "evaluated": evaluated,
        "hits": hits,
        "geneIdColumn": label,
        "hit_ids": hit_ids,
        "evaluable_ids": evaluable,
        "reasons": [],
    }


class CertificateError(ValueError):
    pass


def parse_certificate(data):
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        raise CertificateError("certificate is not valid UTF-8")
    if '"' in text:
        raise CertificateError("certificate must not contain quotes")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if len(lines) != 2:
        raise CertificateError("certificate must have exactly one header and one record")
    if lines[0].split("\t") != QC_COLUMNS:
        raise CertificateError("certificate header columns do not match the contract")
    fields = lines[1].split("\t")
    if len(fields) != len(QC_COLUMNS):
        raise CertificateError("certificate record must have %d fields" % len(QC_COLUMNS))
    (schema, disposition, contrast, gene_col, padj_s, lfc_s, eval_s, hits_s, md5) = fields
    if schema != "1":
        raise CertificateError("unsupported certificate schemaVersion")
    if disposition not in DISPOSITIONS:
        raise CertificateError("unknown disposition")
    if not contrast.strip() or not contrast.isprintable() or len(contrast) > MAX_LABEL:
        raise CertificateError("contrast must be a non-empty bounded single-line string")
    if not gene_col.strip() or not gene_col.isprintable() or len(gene_col) > MAX_LABEL:
        raise CertificateError("geneIdColumn must be a non-empty bounded single-line string")
    if not all(_NUM_RE.fullmatch(n) and n.isascii() for n in (padj_s, lfc_s)):
        raise CertificateError("thresholds are not numeric")
    try:
        padj_t = float(padj_s)
        lfc_t = float(lfc_s)
    except ValueError:
        raise CertificateError("thresholds are not numeric")
    if not (math.isfinite(padj_t) and 0.0 < padj_t <= 1.0):
        raise CertificateError("padjThreshold must satisfy 0 < value <= 1")
    if not (math.isfinite(lfc_t) and lfc_t >= 0.0):
        raise CertificateError("log2fcThreshold must be finite and >= 0")
    if not all(re.fullmatch(r"[0-9]{1,16}", n) for n in (eval_s, hits_s)):
        raise CertificateError("evaluated/hits must be non-negative integers")
    evaluated, hits = int(eval_s), int(hits_s)
    if evaluated > 2**53 - 1 or hits > evaluated:
        raise CertificateError("hits must not exceed evaluated")
    if not _MD5_RE.match(md5):
        raise CertificateError("resultMd5 must be 32 lowercase hex characters")
    expected = "eligible" if evaluated > 0 and hits > 0 else ("skip" if evaluated > 0 else "invalid")
    if disposition != expected or (evaluated == 0 and hits != 0):
        raise CertificateError("disposition is inconsistent with evaluated/hits counts")
    return {
        "disposition": disposition,
        "contrast": contrast,
        "geneIdColumn": gene_col,
        "padjThreshold": padj_t,
        "log2fcThreshold": lfc_t,
        "evaluated": evaluated,
        "hits": hits,
        "resultMd5": md5,
    }


MAX_QC_BYTES = 64 * 1024


def _read_bounded(path, limit):
    with open(path, "rb") as fh:
        data = fh.read(limit + 1)
    if len(data) > limit:
        raise ValueError("file exceeds the %d-byte read limit" % limit)
    return data


def verified_evidence(results_path, qc_path):
    """Read-only verification of a QC certificate against the real results bytes.

    Returns the evaluate_results structure plus contrast/padjThreshold/
    log2fcThreshold/resultMd5 (None when the certificate is unusable).
    Any verification failure yields disposition 'invalid' with reasons.
    """
    cert_fields = {"contrast": None, "padjThreshold": None, "log2fcThreshold": None, "resultMd5": None}

    def fail(reasons):
        out = _invalid(reasons)
        out.update(cert_fields)
        return out

    try:
        qc_data = _read_bounded(qc_path, MAX_QC_BYTES)
    except OSError:
        return fail(["qc certificate not found or unreadable"])
    except ValueError as exc:
        return fail([str(exc)])
    try:
        cert = parse_certificate(qc_data)
    except CertificateError as exc:
        return fail(["malformed qc certificate: %s" % exc])
    cert_fields = {k: cert[k] for k in cert_fields}

    digest = hashlib.md5()
    def hashed_lines(stream):
        for row in stream:
            digest.update(row)
            yield row

    gene_col = None if cert["geneIdColumn"] == UNNAMED_FIRST else cert["geneIdColumn"]
    try:
        with open(results_path, "rb") as stream:
            verdict = evaluate_lines(hashed_lines(stream), gene_col, cert["padjThreshold"], cert["log2fcThreshold"])
    except OSError:
        return fail(["results file not found or unreadable"])
    if digest.hexdigest() != cert["resultMd5"]:
        return fail(["results md5 does not match the certificate (hash drift)"])
    if verdict["disposition"] != "invalid":
        if (verdict["evaluated"] != cert["evaluated"]
                or verdict["hits"] != cert["hits"]
                or verdict["disposition"] != cert["disposition"]
                or verdict["geneIdColumn"] != cert["geneIdColumn"]):
            return fail(["certificate counts/disposition do not match the recomputed results"])
    verdict.update(cert_fields)
    return verdict


def _logaddexp(a, b):
    if a == -math.inf:
        return b
    if b == -math.inf:
        return a
    m = max(a, b)
    return m + math.log(math.exp(a - m) + math.exp(b - m))


def _log_choose(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf_ge(k, N, K, n):
    """One-sided P(X >= k) for X ~ Hypergeometric(N population, K success, n draws).

    Log-space pmf at k plus an upward recurrence, stable for large universes.
    """
    lo = max(0, n - (N - K))
    hi = min(K, n)
    if k <= lo:
        return 1.0
    if k > hi:
        return 0.0
    logp = _log_choose(K, k) + _log_choose(N - K, n - k) - _log_choose(N, n)
    total = logp
    cur = logp
    for i in range(k, hi):
        cur += math.log((K - i) * (n - i)) - math.log((i + 1) * (N - K - n + i + 1))
        total = _logaddexp(total, cur)
    return min(1.0, math.exp(total))


def bh_adjust(pvalues):
    m = len(pvalues)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: (pvalues[i], i))
    q = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvalues[i] * m / rank)
        q[i] = min(1.0, prev)
    return q


class OutputConflict(Exception):
    pass


def write_text_atomic(path, text):
    """Publish a complete table without replacing results from a competing run."""
    data = text.encode("utf-8")
    directory = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(dir=directory, prefix=".deseq2-qc-", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        try:
            os.link(tmp, path)
        except FileExistsError:
            if os.path.islink(path) or not os.path.isfile(path):
                raise OutputConflict(path)
            with open(path, "rb") as fh:
                if fh.read(len(data) + 1) != data:
                    raise OutputConflict(path)
            return "reused"
        return "written"
    finally:
        os.unlink(tmp)
