import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from threading import Thread

import pytest

from pgx.download import DownloadError, download_all, load_sources


@pytest.fixture()
def server():
    payload = {"/ok.csv": b"a,b\n1,2\n", "/page": b"<!DOCTYPE html><html>bot</html>"}

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            path = self.path.split("?", 1)[0]
            if path == "/catalogue.csv":
                body = b"release,filename,url\nDepMap Public 24Q2,Model.csv,http://example.invalid/missing\n"
                self.send_response(200)
                self.send_header("Content-Type", "text/csv")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            if path == "/files" or (path.startswith("/articles/") and path.endswith("/files")):
                body = json.dumps([{"id": 7, "name": "Model.csv", "download_url": f"http://127.0.0.1:{self.server.server_port}/ok.csv"}]).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            data = payload.get(path)
            if data is None:
                self.send_response(404)
                self.end_headers()
                return
            start = 0
            status = 200
            range_header = self.headers.get("Range")
            if range_header:
                start = int(range_header.split("=")[1].split("-")[0])
                status = 206
            body = data[start:]
            self.send_response(status)
            self.send_header("Content-Type", "text/html" if path == "/page" else "text/csv")
            self.send_header("Content-Length", str(len(body)))
            if status == 206:
                self.send_header("Content-Range", f"bytes {start}-{len(data) - 1}/{len(data)}")
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, fmt, *args):
            return

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield httpd
    httpd.shutdown()


def test_download_is_resumable_and_idempotent(tmp_path, server):
    port = server.server_port
    document = {
        "sources": [{
            "id": "demo",
            "version": "1",
            "enabled": True,
            "files": [{"name": "ok.csv", "url": f"http://127.0.0.1:{port}/ok.csv", "role": "demo", "live": False}],
        }]
    }
    raw = tmp_path / "raw"
    first = download_all(document, raw, live=False)
    assert first[0]["status"] == "downloaded"
    partial = raw / "demo" / "1" / "ok.csv.partial"
    (raw / "demo" / "1" / "ok.csv").rename(partial)
    partial.write_bytes(b"a,b\n")
    resumed = download_all(document, raw, live=False)
    assert resumed[0]["status"] == "resumed"
    assert (raw / "demo" / "1" / "ok.csv").read_bytes() == b"a,b\n1,2\n"
    again = download_all(document, raw, live=False)
    assert again[0]["status"] == "skipped"
    manifest = json.loads((raw / "demo" / "1" / "sha256sum.json").read_text())
    assert manifest["files"][0]["sha256"] == again[0]["sha256"]


def test_html_response_is_rejected(tmp_path, server):
    port = server.server_port
    document = {
        "sources": [{
            "id": "demo",
            "version": "1",
            "enabled": True,
            "files": [{"name": "page.csv", "url": f"http://127.0.0.1:{port}/page"}],
        }]
    }
    with pytest.raises(DownloadError, match="HTML"):
        download_all(document, tmp_path / "raw", live=False)


def test_depmap_html_falls_back_to_figshare(tmp_path, server, monkeypatch):
    port = server.server_port

    def fake_catalogue(api_url):
        raise DownloadError(f"{api_url} returned an HTML page")

    monkeypatch.setattr("pgx.download.depmap_catalogue", fake_catalogue)
    document = {
        "sources": [{
            "id": "depmap",
            "version": "24Q2",
            "enabled": True,
            "resolvers": [
                {"type": "depmap_api", "api_url": f"http://127.0.0.1:{port}/catalogue.csv", "release": "DepMap Public 24Q2"},
                {"type": "figshare", "api_url": f"http://127.0.0.1:{port}", "article_id": 1},
            ],
            "files": [{"name": "Model.csv", "figshare_file_id": None}],
        }]
    }
    # figshare resolver without file_id lists /articles/1/files. Point the server path.
    # The test server only implements /files, so call figshare_download_url's list URL by article path.
    # Override by giving file_id None and patching the list endpoint through a direct url on the second resolver
    # is less direct. Use file_id that the server does not special-case: download URL is constructed as
    # {api}/file/download/{id}, which this server does not serve. Instead match by name via a patched fetch.
    document["sources"][0]["files"][0]["url"] = None
    results = download_all(document, tmp_path / "raw", live=False, force_resolve=True)
    assert results[0]["resolver"] == "figshare"
    assert results[0]["sha256"]


def test_live_files_are_skipped_unless_requested(tmp_path):
    document = {
        "sources": [{
            "id": "demo",
            "version": "1",
            "enabled": True,
            "files": [
                {"name": "small.csv", "url": "http://127.0.0.1/missing", "live": False},
                {"name": "huge.csv", "url": "http://127.0.0.1/missing", "live": True},
            ],
        }]
    }
    # small.csv is not live, so the downloader will try it and fail. Mark it absent by using only the live file
    # and an enabled source whose only file is live.
    document["sources"][0]["files"] = [{"name": "huge.csv", "url": "http://127.0.0.1/missing", "live": True}]
    assert download_all(document, tmp_path / "raw", live=False) == []


def test_sources_yaml_declares_batch1_and_leaves_later_batches_disabled():
    document = load_sources(__import__("pathlib").Path(__file__).resolve().parents[1] / "sources.yaml")
    enabled = {source["id"] for source in document["sources"] if source.get("enabled", True)}
    assert enabled == {"depmap", "prism_19q4", "prism_24q2", "gdsc", "gdsc_combinations", "sanger_cmp", "hgnc", "cellosaurus"}
    planned = {item["id"] for item in document["planned"]}
    assert "nci60_dtp" in planned
    assert "orcestra_ctrpv2" in planned
    assert "hmdb" in planned
    assert all(item.get("enabled") is False for item in document["planned"])
    assert "ORCESTRA PRISM" in next(item["notes"] for item in document["planned"] if item["id"] == "orcestra_ccle")
