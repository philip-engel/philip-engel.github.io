"""Small HTTP boundary for the Threefold Explorer Sage engine."""
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import json
import os
import sys
import threading
import subprocess
import traceback

ROOT = Path(__file__).resolve().parent
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


class SageWorker:
    """One warm worker; HTTP health checks never acquire its request lock."""

    def __init__(self):
        self.process = subprocess.Popen(
            [sys.executable, str(ROOT / "sage_worker.py")],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
        self.lock = threading.Lock()
        self.info = json.loads(self.process.stdout.readline())

    def request(self, path, body):
        with self.lock:
            if self.process.poll() is not None:
                raise RuntimeError("The Sage worker stopped.")
            self.process.stdin.write(json.dumps(dict(path=path, body=body))+"\n")
            self.process.stdin.flush()
            line = self.process.stdout.readline()
            if not line:
                raise RuntimeError("The Sage worker stopped during the computation.")
            return json.loads(line)

    def close(self):
        self.process.terminate()
        try:
            self.process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self.process.kill()
            self.process.wait()


WORKER = None


class Handler(BaseHTTPRequestHandler):
    server_version = "ThreefoldExplorer/2"

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
            if WORKER.process.poll() is None:
                self._json(200, WORKER.info)
            else:
                self._json(503, {"status": "unavailable", "error": "The Sage worker stopped."})
        elif self.path == "/api":
            self._json(200, {"name": "Threefold Explorer API", "version": WORKER.info["api_version"]})
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
            response = WORKER.request(self.path, body)
            self._json(response["status"], response["body"])
        except KeyError as error:
            self._json(404 if str(error) == "'Unknown endpoint'" else 422,
                       {"error": str(error).strip("'")})
        except (ValueError, TypeError, ArithmeticError, NotImplementedError) as error:
            self._json(422, {"error": str(error)})
        except (BrokenPipeError, ConnectionResetError):
            self.log_message("Client disconnected before receiving the result.")
        except RuntimeError as error:
            self._json(503, {"error": str(error)})
        except Exception:
            traceback.print_exc()
            self._json(500, {"error": "The Sage computation failed unexpectedly."})

    def log_message(self, format, *args):
        sys.stderr.write("%s - %s\n" % (self.address_string(), format % args))


if __name__ == "__main__":
    WORKER = SageWorker()
    print("Threefold Explorer API on port %d with %d models" %
          (PORT, WORKER.info["models"]), flush=True)
    try:
        ThreadingHTTPServer(("0.0.0.0", PORT), Handler).serve_forever()
    finally:
        WORKER.close()
