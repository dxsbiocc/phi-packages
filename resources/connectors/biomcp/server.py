#!/usr/bin/env python3
"""Run the pinned official BioMCP executable without modifying this package."""

import argparse
import hashlib
import json
import os
import platform
import re
import stat
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

PACKAGE_DIR = Path(__file__).resolve().parent
CHUNK_SIZE = 1024 * 1024


class BootstrapError(Exception):
    """A bootstrap failure that can be reported without a traceback."""


def platform_key(system=None, machine=None):
    system = (system or platform.system()).lower()
    machine = (machine or platform.machine()).lower()
    os_name = {"darwin": "darwin", "macos": "darwin", "linux": "linux"}.get(system)
    arch = {"arm64": "arm64", "aarch64": "arm64", "x86_64": "x64",
            "amd64": "x64", "x64": "x64"}.get(machine)
    key = f"{os_name}-{arch}"
    if key not in {"darwin-arm64", "darwin-x64", "linux-x64"}:
        raise BootstrapError(
            f"unsupported platform {system}/{machine}; supported: macOS arm64/x64, Linux x64"
        )
    return key


def official_url(url):
    parsed = urlsplit(url)
    if (parsed.scheme != "https" or parsed.hostname != "files.pythonhosted.org"
            or parsed.port not in (None, 443) or parsed.username or parsed.password
            or parsed.query or parsed.fragment):
        raise BootstrapError("wheel URL must use HTTPS on files.pythonhosted.org")


def load_upstream(key):
    try:
        metadata = json.loads((PACKAGE_DIR / "upstream.json").read_text(encoding="utf-8"))
        version, pin = metadata["version"], metadata["platforms"][key]
        if metadata["schemaVersion"] != 1 or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", version):
            raise ValueError("invalid version or schema")
        if type(pin["size"]) is not int or pin["size"] <= 0:
            raise ValueError("invalid wheel size")
        if not re.fullmatch(r"[0-9a-f]{64}", pin["sha256"]):
            raise ValueError("invalid wheel SHA256")
        official_url(pin["url"])
        return version, pin
    except (KeyError, TypeError, ValueError) as exc:
        raise BootstrapError(f"invalid upstream.json: {exc}") from exc


def default_cache_dir():
    if platform.system() == "Darwin":
        return Path.home() / "Library" / "Caches" / "Phi" / "connectors" / "biomcp"
    return Path.home() / ".cache" / "Phi" / "connectors" / "biomcp"


def digest_stream(stream):
    digest = hashlib.sha256()
    for chunk in iter(lambda: stream.read(CHUNK_SIZE), b""):
        digest.update(chunk)
    return digest.hexdigest()


def matches(path, size, digest):
    try:
        if path.is_symlink() or not stat.S_ISREG(path.stat().st_mode) or path.stat().st_size != size:
            return False
        with path.open("rb") as stream:
            return digest_stream(stream) == digest
    except OSError:
        return False


def download_wheel(path, pin):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".wheel-", delete=False) as output:
            temporary = Path(output.name)
            request = urllib.request.Request(pin["url"], headers={"User-Agent": "Phi-BioMCP-launcher"})
            with urllib.request.urlopen(request, timeout=30) as response:
                official_url(response.geturl())
                received = 0
                for chunk in iter(lambda: response.read(CHUNK_SIZE), b""):
                    received += len(chunk)
                    if received > pin["size"]:
                        raise BootstrapError("downloaded wheel size exceeds the pinned size")
                    output.write(chunk)
            output.flush()
            os.fsync(output.fileno())
        if temporary.stat().st_size != pin["size"]:
            raise BootstrapError("downloaded wheel size does not match the pinned size")
        if not matches(temporary, pin["size"], pin["sha256"]):
            raise BootstrapError("downloaded wheel SHA256 does not match the pinned checksum")
        os.replace(temporary, path)
    except OSError as exc:
        raise BootstrapError(f"cannot download the pinned BioMCP wheel: {exc}") from exc
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def executable_member(wheel):
    candidates = [info for info in wheel.infolist()
                  if PurePosixPath(info.filename.replace("\\", "/")).name == "biomcp"]
    if len(candidates) != 1:
        raise BootstrapError("wheel must contain exactly one biomcp executable")
    info = candidates[0]
    name = info.orig_filename
    if ("\\" in name or "\x00" in name or ":" in name
            or any(part in ("", ".", "..") for part in name.split("/"))):
        raise BootstrapError("wheel biomcp executable has an unsafe path")
    if info.is_dir() or not stat.S_ISREG(info.external_attr >> 16):
        raise BootstrapError("wheel biomcp executable must be a regular file, not a link")
    return info


def install_executable(wheel, info, path):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".binary-", delete=False) as output:
            temporary = Path(output.name)
            with wheel.open(info) as source:
                for chunk in iter(lambda: source.read(CHUNK_SIZE), b""):
                    output.write(chunk)
            output.flush()
            os.fsync(output.fileno())
            os.fchmod(output.fileno(), 0o700)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def prepare(cache_dir=None):
    key = platform_key()
    version, pin = load_upstream(key)
    root = Path(cache_dir or default_cache_dir()).expanduser().resolve()
    if root == PACKAGE_DIR or PACKAGE_DIR in root.parents:
        raise BootstrapError("cache directory must be outside the connector package")
    cache = root / version / key
    for directory in (root, root / version, cache):
        if directory.is_symlink():
            raise BootstrapError("cache directories must not be symbolic links")
        directory.mkdir(parents=True, exist_ok=True, mode=0o700)
        directory.chmod(0o700)
    import fcntl

    flags = os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0)
    with os.fdopen(os.open(cache / ".lock", flags, 0o600), "a+b") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        archive, binary = cache / "upstream.whl", cache / "biomcp"
        if not matches(archive, pin["size"], pin["sha256"]):
            download_wheel(archive, pin)
        try:
            with zipfile.ZipFile(archive) as wheel:
                info = executable_member(wheel)
                with wheel.open(info) as source:
                    digest = digest_stream(source)
                if not matches(binary, info.file_size, digest):
                    install_executable(wheel, info, binary)
                binary.chmod(0o700)
        except (zipfile.BadZipFile, RuntimeError, NotImplementedError) as exc:
            raise BootstrapError(f"invalid BioMCP wheel: {exc}") from exc
    return binary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true", help="verify and warm the executable cache")
    parser.add_argument("--cache-dir", type=Path, help="override the writable user cache directory")
    args = parser.parse_args(argv)
    try:
        binary = prepare(args.cache_dir)
        if args.prepare:
            print(f"BioMCP cache ready: {binary}", file=sys.stderr)
        else:
            os.execv(str(binary), [str(binary), "serve"])
        return 0
    except (BootstrapError, OSError) as exc:
        print(f"BioMCP launcher: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("BioMCP launcher: interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
