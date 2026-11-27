"""Shared HTTP helpers for omics ingest scanners.

Problem we hit: through the box proxy (DNS → 198.18.0.1), TCP connections can
sit ESTABLISHED with Send-Q stuck. requests' read timeout does not always abort
a blocked send(), so workers freeze for 1–2h despite timeout=60/90.

Fixes:
- socket.setdefaulttimeout covers send/recv
- explicit (connect, read) timeouts
- Connection: close avoids sticky half-dead keep-alives
"""
from __future__ import annotations

import socket
from typing import Any

import requests

# Hard ceiling for any socket op (including write through a wedged proxy).
SOCKET_DEFAULT_TIMEOUT = 60
socket.setdefaulttimeout(SOCKET_DEFAULT_TIMEOUT)

CONNECT_TIMEOUT = 10
READ_TIMEOUT = 45
TIMEOUT = (CONNECT_TIMEOUT, READ_TIMEOUT)

DEFAULT_HEADERS = {
    "Connection": "close",
    "User-Agent": "omics-sniffer/1.0",
}


def http_get(url: str, params: dict | None = None, timeout: Any = TIMEOUT, **kwargs) -> requests.Response:
    headers = {**DEFAULT_HEADERS, **(kwargs.pop("headers", None) or {})}
    return requests.get(url, params=params, timeout=timeout, headers=headers, **kwargs)


def http_post(url: str, data: Any = None, timeout: Any = TIMEOUT, **kwargs) -> requests.Response:
    headers = {**DEFAULT_HEADERS, **(kwargs.pop("headers", None) or {})}
    return requests.post(url, data=data, timeout=timeout, headers=headers, **kwargs)
