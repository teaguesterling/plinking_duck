#!/usr/bin/env python3
"""Static file server with HEAD + single-range GET, logging one TSV line per request.
Usage: range_server.py <root> <port> <logfile>"""
import os, sys, threading
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler

ROOT, PORT, LOG = os.path.abspath(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
lock = threading.Lock()

def log(method, path, status, nbytes):
    with lock, open(LOG, "a") as f:
        f.write(f"{method}\t{path}\t{status}\t{nbytes}\n")

class H(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    def log_message(self, *a):
        pass
    def _file(self):
        p = os.path.normpath(os.path.join(ROOT, self.path.split("?")[0].lstrip("/")))
        return p if p.startswith(ROOT) and os.path.isfile(p) else None
    def _hdr(self, code, length, extra=()):
        self.send_response(code)
        self.send_header("Content-Length", str(length))
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Last-Modified", "Mon, 09 Nov 2026 00:00:00 GMT")
        for k, v in extra:
            self.send_header(k, v)
        self.end_headers()
    def do_HEAD(self):
        p = self._file()
        if not p:
            log("HEAD", self.path, 404, 0); self._hdr(404, 0); return
        log("HEAD", self.path, 200, 0); self._hdr(200, os.path.getsize(p))
    def do_GET(self):
        p = self._file()
        if not p:
            log("GET", self.path, 404, 0); self._hdr(404, 0); return
        size = os.path.getsize(p)
        rng = self.headers.get("Range")
        lo, hi, code = 0, size - 1, 200
        if rng and rng.startswith("bytes="):
            a, b = rng[6:].split(",")[0].split("-")
            lo = int(a) if a else max(0, size - int(b))
            hi = min(int(b), size - 1) if (a and b) else size - 1
            code = 206
        n = max(0, hi - lo + 1)
        extra = [("Content-Range", f"bytes {lo}-{hi}/{size}")] if code == 206 else []
        log("GET", self.path, code, n)
        self._hdr(code, n, extra)
        with open(p, "rb") as f:
            f.seek(lo); self.wfile.write(f.read(n))

ThreadingHTTPServer(("127.0.0.1", PORT), H).serve_forever()
