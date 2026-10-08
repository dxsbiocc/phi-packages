"""Offline security and lifecycle checks for the BioMCP connector launcher."""

import contextlib
import hashlib
import importlib.util
import io
import json
import os
import stat
import tempfile
import time
import unittest
import urllib.error
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch

LAUNCHER = Path(__file__).resolve().parents[1] / "resources" / "connectors" / "biomcp" / "server.py"
spec = importlib.util.spec_from_file_location("biomcp_launcher", LAUNCHER)
launcher = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launcher)

URL = "https://files.pythonhosted.org/packages/fixture/biomcp_cli-0.9.1.whl"
NATIVE = b"offline native executable fixture\n"


def wheel_bytes(members=None):
    members = members if members is not None else [("biomcp_cli.data/scripts/biomcp", stat.S_IFREG | 0o755)]
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w") as wheel:
        for name, mode in members:
            info = zipfile.ZipInfo(name)
            info.create_system = 3
            info.external_attr = mode << 16
            wheel.writestr(info, NATIVE)
    return output.getvalue()


def pin_for(data):
    return {"url": URL, "size": len(data), "sha256": hashlib.sha256(data).hexdigest()}


class Response(io.BytesIO):
    def geturl(self):
        return URL


class PlatformTests(unittest.TestCase):
    def test_supported_platform_aliases(self):
        for system, machine, expected in [
            ("Darwin", "arm64", "darwin-arm64"), ("macOS", "aarch64", "darwin-arm64"),
            ("DARWIN", "AMD64", "darwin-x64"), ("Darwin", "x64", "darwin-x64"),
            ("Darwin", "x86_64", "darwin-x64"), ("Linux", "x86_64", "linux-x64"),
            ("linux", "amd64", "linux-x64"), ("Linux", "x64", "linux-x64"),
        ]:
            with self.subTest(system=system, machine=machine):
                self.assertEqual(launcher.platform_key(system, machine), expected)

    def test_unsupported_platforms(self):
        for system, machine in [("Windows", "amd64"), ("Linux", "aarch64"),
                                ("Linux", "arm64"), ("Darwin", "i386"), ("FreeBSD", "x64")]:
            with self.subTest(system=system, machine=machine):
                with self.assertRaisesRegex(launcher.BootstrapError, "unsupported platform"):
                    launcher.platform_key(system, machine)

    def test_linux_cache_is_stable_across_managed_environments(self):
        expected = Path("/fixture/home/.cache/Phi/connectors/biomcp")
        for environment in ({}, {"XDG_CACHE_HOME": ""},
                            {"XDG_CACHE_HOME": "/runtime/cache/phi-python-1/xdg"}):
            with self.subTest(environment=environment):
                with patch.object(launcher.platform, "system", return_value="Linux"), \
                        patch.object(launcher.Path, "home", return_value=Path("/fixture/home")), \
                        patch.object(launcher.os, "environ", environment):
                    self.assertEqual(launcher.default_cache_dir(), expected)

    def test_official_https_host_is_required(self):
        for url in ["http://files.pythonhosted.org/a.whl", "https://example.com/a.whl",
                    "https://files.pythonhosted.org.evil/a.whl", "https://user@files.pythonhosted.org/a.whl"]:
            with self.subTest(url=url):
                with self.assertRaises(launcher.BootstrapError):
                    launcher.official_url(url)

    def test_packaged_metadata_is_pinned_for_supported_platforms(self):
        for key in ("darwin-arm64", "darwin-x64", "linux-x64"):
            version, pin = launcher.load_upstream(key)
            self.assertEqual(version, "0.9.1")
            self.assertGreater(pin["size"], 0)
            self.assertEqual(len(pin["sha256"]), 64)
            launcher.official_url(pin["url"])

    def test_invalid_metadata_has_clear_error(self):
        invalid = {"schemaVersion": 1, "version": "../0.9.1", "platforms": {"darwin-arm64": pin_for(wheel_bytes())}}
        with patch.object(Path, "read_text", return_value=json.dumps(invalid)):
            with self.assertRaisesRegex(launcher.BootstrapError, "invalid upstream.json"):
                launcher.load_upstream("darwin-arm64")


class CacheTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name).resolve() / "cache"
        self.data = wheel_bytes()
        self.pin = pin_for(self.data)
        self.patches = [patch.object(launcher, "platform_key", return_value="darwin-arm64"),
                        patch.object(launcher, "load_upstream", side_effect=lambda key: ("0.9.1", self.pin)),
                        patch.object(launcher, "default_cache_dir", return_value=self.root),
                        patch.object(launcher.urllib.request, "urlopen", side_effect=AssertionError("unexpected network access"))]
        for item in self.patches:
            item.start()
            self.addCleanup(item.stop)

    @property
    def cache(self):
        return self.root / "0.9.1" / "darwin-arm64"

    def download(self):
        return patch.object(launcher.urllib.request, "urlopen", side_effect=lambda *args, **kwargs: Response(self.data))

    def assert_no_temporary_files(self):
        self.assertFalse(list(self.root.rglob(".wheel-*")))
        self.assertFalse(list(self.root.rglob(".binary-*")))

    def test_verified_download_and_offline_reuse(self):
        with self.download() as download:
            binary = launcher.prepare(self.root)
        download.assert_called_once()
        self.assertEqual(binary.read_bytes(), NATIVE)
        self.assertEqual((self.cache / "upstream.whl").read_bytes(), self.data)
        self.assertEqual(stat.S_IMODE(binary.stat().st_mode), 0o700)
        self.assertEqual(stat.S_IMODE(self.cache.stat().st_mode), 0o700)
        with patch.object(launcher.urllib.request, "urlopen", side_effect=AssertionError("must stay offline")):
            self.assertEqual(launcher.prepare(self.root), binary)
        self.assert_no_temporary_files()

    def test_corrupted_executable_is_repaired_from_verified_wheel(self):
        with self.download():
            binary = launcher.prepare(self.root)
        binary.write_bytes(b"untrusted replacement")
        binary.chmod(0o777)
        with patch.object(launcher.urllib.request, "urlopen", side_effect=AssertionError("must stay offline")):
            launcher.prepare(self.root)
        self.assertEqual(binary.read_bytes(), NATIVE)
        self.assertEqual(stat.S_IMODE(binary.stat().st_mode), 0o700)
        self.assert_no_temporary_files()

    def test_executable_symlink_is_replaced_without_touching_target(self):
        with self.download():
            binary = launcher.prepare(self.root)
        target = Path(self.directory.name) / "target"
        target.write_bytes(NATIVE)
        binary.unlink()
        binary.symlink_to(target)
        with patch.object(launcher.urllib.request, "urlopen", side_effect=AssertionError("must stay offline")):
            launcher.prepare(self.root)
        self.assertFalse(binary.is_symlink())
        self.assertEqual(target.read_bytes(), NATIVE)
        self.assertEqual(binary.read_bytes(), NATIVE)

    def test_corrupted_wheel_is_replaced_by_verified_download(self):
        with self.download():
            launcher.prepare(self.root)
        (self.cache / "upstream.whl").write_bytes(b"untrusted wheel")
        with self.download() as download:
            launcher.prepare(self.root)
        download.assert_called_once()
        self.assertEqual((self.cache / "upstream.whl").read_bytes(), self.data)

    def test_corrupted_wheel_cannot_authorize_cached_executable_offline(self):
        with self.download():
            launcher.prepare(self.root)
        (self.cache / "upstream.whl").write_bytes(b"untrusted wheel")
        with patch.object(launcher.urllib.request, "urlopen", side_effect=urllib.error.URLError("offline")):
            with self.assertRaisesRegex(launcher.BootstrapError, "cannot download"):
                launcher.prepare(self.root)
        self.assert_no_temporary_files()

    def test_checksum_failure_occurs_before_zip_contents_are_used(self):
        self.pin["sha256"] = "0" * 64
        with self.download(), patch.object(launcher.zipfile, "ZipFile") as archive:
            with self.assertRaisesRegex(launcher.BootstrapError, "SHA256"):
                launcher.prepare(self.root)
        archive.assert_not_called()
        self.assertFalse((self.cache / "upstream.whl").exists())
        self.assert_no_temporary_files()

    def test_short_and_oversized_downloads_are_rejected(self):
        for delta in (-1, 1):
            with self.subTest(delta=delta):
                self.pin["size"] = len(self.data) + delta
                with self.download():
                    with self.assertRaisesRegex(launcher.BootstrapError, "size"):
                        launcher.prepare(self.root)
                self.assertFalse((self.cache / "upstream.whl").exists())
                self.assert_no_temporary_files()

    def test_partial_download_failure_preserves_existing_archive(self):
        self.cache.mkdir(parents=True)
        archive = self.cache / "upstream.whl"
        archive.write_bytes(self.data)

        class BrokenResponse(Response):
            def read(self, size):
                if self.tell():
                    raise OSError("connection dropped")
                return super().read(10)

        with patch.object(launcher.urllib.request, "urlopen", return_value=BrokenResponse(self.data)):
            with self.assertRaisesRegex(launcher.BootstrapError, "connection dropped"):
                launcher.download_wheel(archive, self.pin)
        self.assertEqual(archive.read_bytes(), self.data)
        self.assert_no_temporary_files()

    def test_redirect_to_unapproved_host_is_rejected(self):
        response = Response(self.data)
        response.geturl = lambda: "https://example.com/wheel.whl"
        with patch.object(launcher.urllib.request, "urlopen", return_value=response):
            with self.assertRaisesRegex(launcher.BootstrapError, "HTTPS on files.pythonhosted.org"):
                launcher.prepare(self.root)
        self.assert_no_temporary_files()

    def test_cache_inside_package_is_rejected(self):
        with patch.object(launcher, "PACKAGE_DIR", Path(self.directory.name).resolve()):
            with self.assertRaisesRegex(launcher.BootstrapError, "outside the connector package"):
                launcher.prepare(self.root)
        self.assertFalse(self.root.exists())

    def test_concurrent_starts_download_once(self):
        def response(*args, **kwargs):
            time.sleep(0.05)
            return Response(self.data)

        with patch.object(launcher.urllib.request, "urlopen", side_effect=response) as download:
            with ThreadPoolExecutor(max_workers=2) as pool:
                binaries = list(pool.map(lambda _: launcher.prepare(self.root), range(2)))
        self.assertEqual(binaries[0], binaries[1])
        download.assert_called_once()
        self.assert_no_temporary_files()

    def test_unsafe_missing_and_ambiguous_executable_members(self):
        regular = stat.S_IFREG | 0o755
        fixtures = [
            [], [("other", regular)], [("../biomcp", regular)], [("/biomcp", regular)],
            [("bin\\biomcp", regular)], [("bin//biomcp", regular)], [("bin/./biomcp", regular)],
            [("C:/biomcp", regular)], [("bin/biomcp/", stat.S_IFDIR | 0o755)],
            [("bin/biomcp", stat.S_IFLNK | 0o777)], [("bin/biomcp", stat.S_IFIFO | 0o700)],
            [("a/biomcp", regular), ("b/biomcp", regular)],
        ]
        for members in fixtures:
            with self.subTest(members=members):
                self.data = wheel_bytes(members)
                self.pin = pin_for(self.data)
                with self.download():
                    with self.assertRaises(launcher.BootstrapError):
                        launcher.prepare(self.root)
                self.assertFalse((self.cache / "biomcp").exists())
                self.assert_no_temporary_files()

    def test_only_executable_is_extracted(self):
        self.data = wheel_bytes([("scripts/biomcp", stat.S_IFREG | 0o755),
                                 ("../do-not-extract", stat.S_IFREG | 0o644)])
        self.pin = pin_for(self.data)
        with self.download():
            launcher.prepare(self.root)
        self.assertEqual({item.name for item in self.cache.iterdir()}, {".lock", "upstream.whl", "biomcp"})
        self.assertFalse((self.root / "do-not-extract").exists())

    def test_default_launch_execs_only_serve_without_stdout(self):
        stdout, stderr = io.StringIO(), io.StringIO()
        with self.download(), patch.object(launcher.os, "execv") as execute:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                self.assertEqual(launcher.main([]), 0)
        binary = str(self.cache / "biomcp")
        execute.assert_called_once_with(binary, [binary, "serve"])
        self.assertEqual(stdout.getvalue(), "")
        self.assertEqual(stderr.getvalue(), "")

    def test_prepare_flag_reports_only_on_stderr(self):
        stdout, stderr = io.StringIO(), io.StringIO()
        with self.download(), patch.object(launcher.os, "execv") as execute:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                self.assertEqual(launcher.main(["--prepare", "--cache-dir", str(self.root)]), 0)
        execute.assert_not_called()
        self.assertEqual(stdout.getvalue(), "")
        self.assertIn("BioMCP cache ready:", stderr.getvalue())

    def test_bootstrap_error_reports_only_on_stderr(self):
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch.object(launcher.urllib.request, "urlopen", side_effect=urllib.error.URLError("offline")):
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                self.assertEqual(launcher.main([]), 1)
        self.assertEqual(stdout.getvalue(), "")
        self.assertIn("cannot download", stderr.getvalue())
        self.assert_no_temporary_files()

    def test_extra_native_arguments_are_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()), patch.object(launcher, "prepare") as prepare:
            with self.assertRaises(SystemExit) as error:
                launcher.main(["--arbitrary-native-flag"])
        self.assertEqual(error.exception.code, 2)
        prepare.assert_not_called()


if __name__ == "__main__":
    unittest.main()
