#!/usr/bin/env python3
"""Local PDF text extraction helpers (pdftotext + optional OCR).

Stdlib + host tools only. Soft-fail when tools are missing — never hard-fail
the whole ingest solely because OCR/TextIn is absent. TextIn is NOT a required
parser and is NOT the KB backbone.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Optional


MIN_USEFUL_CHARS = 40
OCR_DPI = 200
MAX_OCR_PAGES = 40


def which_pdftotext() -> Optional[str]:
    return shutil.which("pdftotext")


def which_tesseract() -> Optional[str]:
    return shutil.which("tesseract")


def extract_with_pdftotext(path: Path, *, pdftotext_bin: Optional[str] = None) -> dict[str, Any]:
    """Run pdftotext -layout -enc UTF-8. Returns {ok, text, tool, warning?}."""
    bin_path = pdftotext_bin or which_pdftotext()
    if not bin_path:
        return {
            "ok": False,
            "text": "",
            "tool": None,
            "warning": "pdftotext_missing",
            "message": "pdftotext not on PATH; install poppler-utils or use a text sidecar",
        }
    try:
        proc = subprocess.run(
            [bin_path, "-layout", "-enc", "UTF-8", str(path), "-"],
            capture_output=True,
            timeout=120,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as e:
        return {
            "ok": False,
            "text": "",
            "tool": bin_path,
            "warning": "pdftotext_failed",
            "message": str(e),
        }
    text = (proc.stdout or b"").decode("utf-8", errors="replace")
    if proc.returncode != 0 and not text.strip():
        err = (proc.stderr or b"").decode("utf-8", errors="replace")[:500]
        return {
            "ok": False,
            "text": "",
            "tool": bin_path,
            "warning": "pdftotext_nonzero",
            "message": err or f"exit {proc.returncode}",
        }
    return {
        "ok": True,
        "text": text,
        "tool": bin_path,
        "n_chars": len(text),
        "useful": len(text.strip()) >= MIN_USEFUL_CHARS,
    }


def extract_with_tesseract(
    path: Path,
    *,
    tesseract_bin: Optional[str] = None,
    lang: str = "eng",
) -> dict[str, Any]:
    """Optional OCR via tesseract (soft-fail if missing)."""
    bin_path = tesseract_bin or which_tesseract()
    if not bin_path:
        return {
            "ok": False,
            "text": "",
            "tool": None,
            "warning": "tesseract_missing",
            "message": (
                "tesseract not on PATH; OCR skipped. Install tesseract-ocr or set "
                "LIT_PDF_OCR_CMD / --ocr_cmd for an external OCR hook."
            ),
        }
    # tesseract reads images, not PDFs: rasterize pages with pdftoppm (poppler) first.
    pdftoppm = shutil.which("pdftoppm")
    if not pdftoppm:
        return {
            "ok": False,
            "text": "",
            "tool": bin_path,
            "warning": "pdftoppm_missing",
            "message": "pdftoppm (poppler) not on PATH; cannot rasterize the PDF for OCR.",
        }
    with tempfile.TemporaryDirectory() as td:
        try:
            proc = subprocess.run(
                [pdftoppm, "-r", str(OCR_DPI), "-l", str(MAX_OCR_PAGES), "-gray", "-png", str(path), str(Path(td) / "p")],
                capture_output=True,
                timeout=180,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as e:
            return {"ok": False, "text": "", "tool": pdftoppm, "warning": "pdftoppm_failed", "message": str(e)}
        pages = sorted(Path(td).glob("p-*.png"))
        if not pages:
            err = (proc.stderr or b"").decode("utf-8", errors="replace")[:300]
            return {"ok": False, "text": "", "tool": pdftoppm, "warning": "pdftoppm_no_pages", "message": err or f"exit {proc.returncode}"}
        chunks: list[str] = []
        for img in pages:
            try:
                r = subprocess.run(
                    [bin_path, str(img), "stdout", "-l", lang],
                    capture_output=True,
                    timeout=120,
                    check=False,
                )
            except (OSError, subprocess.TimeoutExpired) as e:
                return {"ok": False, "text": "\n\n".join(chunks), "tool": bin_path, "warning": "tesseract_failed", "message": str(e)}
            if r.returncode == 0:
                chunks.append(r.stdout.decode("utf-8", errors="replace").strip())
        text = "\n\n".join(c for c in chunks if c)
        if not text.strip():
            return {"ok": False, "text": "", "tool": bin_path, "warning": "tesseract_no_output", "message": f"{len(pages)} page(s) rendered, no text recognised"}
        out = {
            "ok": True,
            "text": text,
            "tool": bin_path,
            "n_chars": len(text),
            "n_pages_ocr": len(pages),
            "useful": len(text.strip()) >= MIN_USEFUL_CHARS,
        }
        if len(pages) >= MAX_OCR_PAGES:
            out["message"] = f"OCR limited to the first {MAX_OCR_PAGES} pages"
        return out


def extract_with_ocr_cmd(path: Path, ocr_cmd: str) -> dict[str, Any]:
    """External OCR hook: shell template with {path} and optional {outdir}.

    Example: LIT_PDF_OCR_CMD='my-ocr --in {path} --out {outdir}/out.txt'
    Must write UTF-8 text to {outdir}/out.txt OR print text on stdout.
    """
    with tempfile.TemporaryDirectory() as td:
        outdir = Path(td)
        out_txt = outdir / "out.txt"
        cmd = ocr_cmd.format(path=str(path), outdir=str(outdir), out=str(out_txt))
        try:
            proc = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                timeout=300,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as e:
            return {
                "ok": False,
                "text": "",
                "tool": "ocr_cmd",
                "warning": "ocr_cmd_failed",
                "message": str(e),
                "ocr_cmd": ocr_cmd,
            }
        text = ""
        if out_txt.is_file():
            text = out_txt.read_text(encoding="utf-8", errors="replace")
        elif proc.stdout:
            text = proc.stdout.decode("utf-8", errors="replace")
        if not text.strip():
            err = (proc.stderr or b"").decode("utf-8", errors="replace")[:500]
            return {
                "ok": False,
                "text": "",
                "tool": "ocr_cmd",
                "warning": "ocr_cmd_empty",
                "message": err or f"exit {proc.returncode}",
                "ocr_cmd": ocr_cmd,
            }
        return {
            "ok": True,
            "text": text,
            "tool": "ocr_cmd",
            "n_chars": len(text),
            "useful": len(text.strip()) >= MIN_USEFUL_CHARS,
            "ocr_cmd": ocr_cmd,
        }


def extract_pdf_text(
    path: Path,
    *,
    try_ocr: bool = True,
    ocr_cmd: Optional[str] = None,
    pdftotext_bin: Optional[str] = None,
    tesseract_bin: Optional[str] = None,
) -> dict[str, Any]:
    """Best-effort text from a local PDF. Never raises for missing OCR tools.

    Order: pdftotext → (if empty/near-empty and try_ocr) ocr_cmd or tesseract.
    """
    path = Path(path)
    warnings: list[str] = []
    steps: list[dict[str, Any]] = []

    digital = extract_with_pdftotext(path, pdftotext_bin=pdftotext_bin)
    steps.append({"step": "pdftotext", **{k: digital.get(k) for k in ("ok", "warning", "n_chars", "useful", "message")}})
    text = digital.get("text") or ""
    if digital.get("ok") and digital.get("useful"):
        return {
            "ok": True,
            "text": text,
            "method": "pdftotext",
            "tool": digital.get("tool"),
            "n_chars": len(text),
            "warnings": warnings,
            "steps": steps,
        }
    if digital.get("warning"):
        warnings.append(str(digital["warning"]))

    if not try_ocr:
        return {
            "ok": bool(text.strip()),
            "text": text,
            "method": "pdftotext" if text.strip() else None,
            "tool": digital.get("tool"),
            "n_chars": len(text),
            "warnings": warnings + (["near_empty_no_ocr"] if text.strip() else ["empty_no_ocr"]),
            "steps": steps,
        }

    hook = ocr_cmd or os.environ.get("LIT_PDF_OCR_CMD") or ""
    if hook.strip():
        ocr = extract_with_ocr_cmd(path, hook.strip())
        steps.append({"step": "ocr_cmd", **{k: ocr.get(k) for k in ("ok", "warning", "n_chars", "useful", "message")}})
        if ocr.get("ok") and (ocr.get("text") or "").strip():
            return {
                "ok": True,
                "text": ocr["text"],
                "method": "ocr_cmd",
                "tool": "ocr_cmd",
                "n_chars": len(ocr["text"]),
                "warnings": warnings,
                "steps": steps,
            }
        if ocr.get("warning"):
            warnings.append(str(ocr["warning"]))
    else:
        ocr = extract_with_tesseract(path, tesseract_bin=tesseract_bin)
        steps.append({"step": "tesseract", **{k: ocr.get(k) for k in ("ok", "warning", "n_chars", "useful", "message")}})
        if ocr.get("ok") and (ocr.get("text") or "").strip():
            return {
                "ok": True,
                "text": ocr["text"],
                "method": "tesseract",
                "tool": ocr.get("tool"),
                "n_chars": len(ocr["text"]),
                "warnings": warnings,
                "steps": steps,
            }
        if ocr.get("warning"):
            warnings.append(str(ocr["warning"]))
            if ocr.get("message"):
                warnings.append(f"ocr_message:{ocr['message'][:200]}")

    # Fall back to whatever digital text we got (may be empty)
    return {
        "ok": bool(text.strip()),
        "text": text,
        "method": "pdftotext" if text.strip() else None,
        "tool": digital.get("tool"),
        "n_chars": len(text),
        "warnings": warnings or ["no_text_extracted"],
        "steps": steps,
    }


def first_nonzero_title_line(text: str, *, fallback: str) -> str:
    """Pick a plausible title from early non-empty lines."""
    for line in (text or "").splitlines()[:40]:
        s = line.strip()
        if len(s) < 8:
            continue
        if s.lower().startswith("%pdf"):
            continue
        if s.isdigit():
            continue
        return s[:500]
    return fallback


_JUNK_TITLE = re.compile(r"^(microsoft word|untitled|document\d*|\w+\.(docx?|pdf|tex)$)", re.I)


def pdf_info_title(path: Path, *, pdfinfo_bin: Optional[str] = None) -> Optional[str]:
    """Title from the PDF's own metadata (pdfinfo), or None when absent or not a real title."""
    bin_path = pdfinfo_bin or shutil.which("pdfinfo")
    if not bin_path:
        return None
    try:
        proc = subprocess.run([bin_path, str(path)], capture_output=True, timeout=30, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return None
    for line in proc.stdout.decode("utf-8", errors="replace").splitlines():
        if line.lower().startswith("title:"):
            title = re.sub(r"\s+", " ", line.split(":", 1)[1]).strip()
            if len(title) >= 20 and not _JUNK_TITLE.match(title) and len(title.split()) >= 5:
                return title[:500]
    return None


_ABSTRACT_RE = re.compile(
    r"\babstract\b[\s:.|—–-]*(?P<body>.{200,3500}?)"
    r"(?=\n\s*(?:keywords?|key words|introduction|1[.\s]+introduction|background)\b|\Z)",
    re.I | re.S,
)


def guess_abstract(text: str) -> Optional[str]:
    """The paragraph after an 'Abstract' heading in the first pages, or None."""
    head = (text or "")[:15000]
    m = _ABSTRACT_RE.search(head)
    if not m:
        return None
    body = re.sub(r"\s*\n\s*", " ", m.group("body")).strip()
    return body if len(body) >= 200 else None


def front_matter_text(path: Path, *, pages: int = 2, pdftotext_bin: Optional[str] = None) -> str:
    """Reading-order text of the first pages, for finding the paper's own identifiers.

    Unlike the -layout text used for the body, two-column pages are not interleaved, so
    wrapped DOIs such as "10.1038/" + "s41575-020-0342-4" stay intact. Empty when unavailable.
    """
    bin_path = pdftotext_bin or which_pdftotext()
    if not bin_path:
        return ""
    try:
        proc = subprocess.run(
            [bin_path, "-f", "1", "-l", str(pages), "-enc", "UTF-8", str(path), "-"],
            capture_output=True, timeout=60, check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    return proc.stdout.decode("utf-8", errors="replace")
