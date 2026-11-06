#!/usr/bin/env python3
"""Local PDF text extraction helpers (pdftotext + optional OCR).

Stdlib + host tools only. Soft-fail when tools are missing — never hard-fail
the whole ingest solely because OCR/TextIn is absent. TextIn is NOT a required
parser and is NOT the KB backbone.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Optional


MIN_USEFUL_CHARS = 40


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
    # tesseract can OCR PDF on many builds; if not, soft-fail
    with tempfile.TemporaryDirectory() as td:
        out_base = Path(td) / "ocr"
        try:
            proc = subprocess.run(
                [bin_path, str(path), str(out_base), "-l", lang],
                capture_output=True,
                timeout=180,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as e:
            return {
                "ok": False,
                "text": "",
                "tool": bin_path,
                "warning": "tesseract_failed",
                "message": str(e),
            }
        txt_path = Path(str(out_base) + ".txt")
        if not txt_path.is_file():
            err = (proc.stderr or b"").decode("utf-8", errors="replace")[:500]
            return {
                "ok": False,
                "text": "",
                "tool": bin_path,
                "warning": "tesseract_no_output",
                "message": err or f"exit {proc.returncode}",
            }
        text = txt_path.read_text(encoding="utf-8", errors="replace")
        return {
            "ok": True,
            "text": text,
            "tool": bin_path,
            "n_chars": len(text),
            "useful": len(text.strip()) >= MIN_USEFUL_CHARS,
        }


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
