"""Idempotent, resumable downloads into raw/<source>/<version>/."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import yaml

USER_AGENT = "phi-pgx/1.0 (+https://github.com/dxsbiocc/phi-packages)"
CHUNK = 1024 * 1024


class DownloadError(RuntimeError):
    pass


def load_sources(path: Path) -> dict:
    with path.open() as handle:
        document = yaml.safe_load(handle)
    if not isinstance(document, dict) or "sources" not in document:
        raise DownloadError(f"{path} has no sources list")
    return document


def live_enabled(explicit: bool | None = None) -> bool:
    if explicit is not None:
        return explicit
    return os.environ.get("PHI_PGX_LIVE", "").strip() in {"1", "true", "TRUE", "yes"}


def _request(url: str, headers: dict[str, str] | None = None, timeout: int = 120):
    merged = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    if headers:
        merged.update(headers)
    return urlopen(Request(url, headers=merged), timeout=timeout)


def _looks_like_html(payload: bytes, content_type: str) -> bool:
    if "html" in content_type.lower():
        return True
    sample = payload[:200].lstrip().lower()
    return sample.startswith(b"<!doctype html") or sample.startswith(b"<html")


def fetch_bytes(url: str, timeout: int = 120) -> tuple[bytes, str]:
    with _request(url, timeout=timeout) as response:
        content_type = response.headers.get("Content-Type", "")
        payload = response.read()
    if _looks_like_html(payload, content_type):
        raise DownloadError(f"{url} returned an HTML page, not a data file. DepMap portal downloads are behind a bot check from some networks; set a figshare or direct URL.")
    return payload, content_type


def depmap_catalogue(api_url: str) -> list[dict[str, str]]:
    payload, _content_type = fetch_bytes(api_url)
    text = payload.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    if not reader.fieldnames or "filename" not in reader.fieldnames:
        raise DownloadError(f"{api_url} did not return a filename catalogue")
    return [dict(row) for row in reader]


def figshare_download_url(api_url: str, article_id: int, filename: str, file_id: int | None) -> str:
    if file_id:
        return f"{api_url.rstrip('/')}/file/download/{file_id}"
    page = 1
    while True:
        url = f"{api_url.rstrip('/')}/articles/{article_id}/files?page={page}&page_size=100"
        payload, _content_type = fetch_bytes(url)
        batch = json.loads(payload)
        if not isinstance(batch, list):
            raise DownloadError(f"figshare article {article_id} did not return a file list")
        for item in batch:
            if item.get("name") == filename and item.get("download_url"):
                return str(item["download_url"])
        if len(batch) < 100:
            break
        page += 1
    raise DownloadError(f"{filename} is not in figshare article {article_id}")


def resolve_url(source: dict, spec: dict, force: bool = False) -> tuple[str, str]:
    configured = spec.get("url")
    if configured and not force:
        return str(configured), "configured_url"
    errors: list[str] = []
    for resolver in source.get("resolvers") or []:
        kind = resolver.get("type")
        try:
            if kind == "direct":
                if not resolver.get("url"):
                    raise DownloadError("direct resolver has no url")
                return str(resolver["url"]), "direct"
            if kind == "depmap_api":
                rows = depmap_catalogue(str(resolver["api_url"]))
                release = str(resolver["release"])
                for row in rows:
                    if row.get("release") == release and row.get("filename") == spec["name"] and row.get("url"):
                        return str(row["url"]), "depmap_api"
                raise DownloadError(f"no url for {release} {spec['name']}")
            if kind == "figshare":
                file_id = spec.get("figshare_file_id")
                return figshare_download_url(
                    str(resolver.get("api_url") or "https://api.figshare.com/v2"),
                    int(resolver["article_id"]),
                    str(spec["name"]),
                    int(file_id) if file_id else None,
                ), "figshare"
            raise DownloadError(f"unknown resolver {kind}")
        except (DownloadError, HTTPError, URLError, TimeoutError, json.JSONDecodeError, UnicodeError) as exc:
            errors.append(f"{kind}: {exc}")
    detail = "; ".join(errors) or "no resolvers configured"
    raise DownloadError(f"could not resolve {spec['name']}: {detail}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _manifest_path(directory: Path) -> Path:
    return directory / "sha256sum.json"


def read_manifest(directory: Path) -> dict:
    path = _manifest_path(directory)
    if not path.exists():
        return {"files": []}
    return json.loads(path.read_text())


def write_manifest(directory: Path, source: str, version: str, files: list[dict]) -> None:
    payload = {
        "source": source,
        "version": version,
        "files": files,
    }
    _manifest_path(directory).write_text(json.dumps(payload, indent=2) + "\n")


def download_file(url: str, dest: Path, expected_sha256: str | None = None) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        actual = sha256_file(dest)
        if expected_sha256 is None or actual == expected_sha256:
            return {"name": dest.name, "sha256": actual, "bytes": dest.stat().st_size, "url": url, "status": "skipped", "retrieved_at": datetime.now(timezone.utc).isoformat()}
        raise DownloadError(f"{dest.name} sha256 {actual} does not match {expected_sha256}")
    partial = dest.with_suffix(dest.suffix + ".partial")
    existing = partial.stat().st_size if partial.exists() else 0
    headers = {"Range": f"bytes={existing}-"} if existing else {}
    try:
        response = _request(url, headers=headers)
    except HTTPError as exc:
        if exc.code == 416 and existing:
            partial.replace(dest)
            return _finish(dest, url, expected_sha256, "resumed")
        raise DownloadError(f"{url} returned HTTP {exc.code}") from exc
    with response:
        status = getattr(response, "status", 200)
        content_type = response.headers.get("Content-Type", "")
        if status == 200 and existing:
            partial.unlink(missing_ok=True)
            existing = 0
        mode = "ab" if status == 206 and existing else "wb"
        if mode == "wb" and partial.exists():
            partial.unlink()
        with partial.open(mode) as handle:
            first = response.read(CHUNK)
            if not existing and _looks_like_html(first, content_type):
                partial.unlink(missing_ok=True)
                raise DownloadError(f"{url} returned HTML. Refusing to store a bot-check page as {dest.name}.")
            if first:
                handle.write(first)
            for chunk in iter(lambda: response.read(CHUNK), b""):
                handle.write(chunk)
    partial.replace(dest)
    return _finish(dest, url, expected_sha256, "downloaded" if not existing else "resumed")


def _finish(dest: Path, url: str, expected_sha256: str | None, status: str) -> dict:
    actual = sha256_file(dest)
    if expected_sha256 and actual != expected_sha256:
        dest.unlink(missing_ok=True)
        raise DownloadError(f"{dest.name} sha256 {actual} does not match {expected_sha256}")
    return {
        "name": dest.name,
        "sha256": actual,
        "bytes": dest.stat().st_size,
        "url": url,
        "status": status,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
    }


def download_source(source: dict, raw_root: Path, live: bool, force_resolve: bool = False, only: set[str] | None = None) -> list[dict]:
    if not source.get("enabled", True):
        return []
    if only and source["id"] not in only:
        return []
    version = str(source["version"])
    directory = raw_root / source["id"] / version
    directory.mkdir(parents=True, exist_ok=True)
    previous = {item["name"]: item for item in read_manifest(directory).get("files", [])}
    results: list[dict] = []
    for spec in source.get("files") or []:
        if spec.get("live") and not live:
            continue
        url, resolver = resolve_url(source, spec, force=force_resolve)
        record = download_file(url, directory / spec["name"], spec.get("sha256"))
        record["resolver"] = resolver
        record["role"] = spec.get("role")
        previous[spec["name"]] = record
        results.append(record)
    write_manifest(directory, source["id"], version, list(previous.values()))
    return results


def download_all(document: dict, raw_root: Path, live: bool, force_resolve: bool = False, only: set[str] | None = None) -> list[dict]:
    results: list[dict] = []
    for source in document.get("sources") or []:
        results.extend(download_source(source, raw_root, live, force_resolve, only))
    return results
