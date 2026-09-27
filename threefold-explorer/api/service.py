"""Small HTTP boundary for the Threefold Explorer Sage engine."""
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
import json
import os
import sys
import threading

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "engine"))

from explorer_api import (
    compute,
    describe_os_entry,
    log_transform_schema,
    section_and_linearization_schema,
)
from local_model_database import LocalModelDatabase

PORT = int(os.environ.get("PORT", "8000"))
MAX_BODY = 128 * 1024
DEFAULT_ORIGINS = (
    "https://philip-engel.github.io",
    "http://127.0.0.1:8080",
    "http://localhost:8080",
)
ALLOWED_ORIGINS = {
    value.strip().rstrip("/")
    for value in os.environ.get("ALLOWED_ORIGINS", ",".join(DEFAULT_ORIGINS)).split(",")
    if value.strip()
}
DATABASE = LocalModelDatabase()
COMPUTE_LOCK = threading.Lock()


def payload_for(path, body):
    if path == "/api/os-entry":
        return describe_os_entry(body["os_entry"], profile=body.get("profile", "default"))
    if path == "/api/sections":
        return section_and_linearization_schema(
            body["os_entry"], body["P"], body["Q"],
            profile=body.get("profile", "default"),
            smooth_slots=body.get("smooth_slots", 0),
        )
    if path == "/api/log-schema":
        return log_transform_schema(
            body["os_entry"], body["P"], body["Q"], body["linearization_divisor"],
            profile=body.get("profile", "default"),
        )
    if path == "/api/compute":
        with COMPUTE_LOCK:
            return compute(body, verbose=False, database=DATABASE)
    raise KeyError("Unknown endpoint")


class Handler(BaseHTTPRequestHandler):
    server_version = "ThreefoldExplorer/1"

    def _origin(self):
        origin = self.headers.get("Origin", "").rstrip("/")
        return origin if origin in ALLOWED_ORIGINS else None

    def _headers(self, status, length):
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(length))
        self.send_header("Cache-Control", "no-store")
        origin = self._origin()
        if origin:
            self.send_header("Access-Control-Allow-Origin", origin)
            self.send_header("Vary", "Origin")
        self.end_headers()

    def _json(self, status, value):
        data = json.dumps(value, separators=(",", ":")).encode("utf-8")
        self._headers(status, len(data))
        self.wfile.write(data)

    def do_OPTIONS(self):
        origin = self._origin()
        if not origin:
            self._json(403, {"error": "Origin is not allowed."})
            return
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", origin)
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Max-Age", "600")
        self.send_header("Vary", "Origin")
        self.end_headers()

    def do_GET(self):
        if self.path == "/health":
            self._json(200, {"status": "ok", "models": len(DATABASE.index["models"])})
        elif self.path == "/api":
            self._json(200, {"name": "Threefold Explorer API", "version": 1})
        else:
            self._json(404, {"error": "Not found."})

    def do_POST(self):
        if self._origin() is None and self.headers.get("Origin"):
            self._json(403, {"error": "Origin is not allowed."})
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > MAX_BODY:
                raise ValueError("Request body must be between 1 byte and 128 KB.")
            body = json.loads(self.rfile.read(length))
            if not isinstance(body, dict):
                raise ValueError("The request body must be a JSON object.")
            self._json(200, payload_for(self.path, body))
        except KeyError as error:
            self._json(404 if str(error) == "'Unknown endpoint'" else 422,
                       {"error": str(error).strip("'")})
        except (ValueError, TypeError, ArithmeticError, NotImplementedError) as error:
            self._json(422, {"error": str(error)})
        except Exception:
            self._json(500, {"error": "The Sage computation failed unexpectedly."})

    def log_message(self, format, *args):
        sys.stderr.write("%s - %s\n" % (self.address_string(), format % args))


if __name__ == "__main__":
    print("Threefold Explorer API on port %d with %d models" %
          (PORT, len(DATABASE.index["models"])), flush=True)
    HTTPServer(("0.0.0.0", PORT), Handler).serve_forever()
