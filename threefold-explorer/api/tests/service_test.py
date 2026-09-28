"""Integration check: health remains responsive during a real Sage computation.

Run with Sage's Python from any directory. Starts its own temporary HTTP service.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import os
import socket
import subprocess
import sys
import tempfile
import time
import urllib.request
import urllib.error

API = Path(__file__).resolve().parents[1]
with socket.socket() as sock:
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
base = "http://127.0.0.1:%d" % port

def request(path, body=None):
    data = None if body is None else json.dumps(body).encode()
    req = urllib.request.Request(base+path, data=data,
                                 headers={"Content-Type":"application/json"})
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.load(response)

with tempfile.TemporaryFile(mode="w+") as logs:
    env = dict(os.environ, PORT=str(port))
    server = subprocess.Popen([sys.executable, str(API/"service.py")], env=env,
                              stdout=logs, stderr=logs)
    try:
        deadline = time.monotonic()+45
        while True:
            try:
                assert request("/health")["models"] == 865
                break
            except (OSError, AssertionError):
                if server.poll() is not None or time.monotonic()>deadline:
                    raise
                time.sleep(.2)
        with ThreadPoolExecutor(max_workers=1) as pool:
            job = pool.submit(request, "/api/compute",
                dict(os_entry=56,P=[30],Q=[1],linearization_divisor=[0,0,0,0,1,0],
                     log_data=[None]*5+[[0,0,0,1]],coordinates="invariant"))
            delays = []
            while not job.done():
                started = time.monotonic()
                assert request("/health")["api_version"] == 2
                delays.append(time.monotonic()-started)
                time.sleep(.1)
            result = job.result()
        assert result["S6_for_supplied_smooth_model"]
        assert result["fundamental_group"]["trivial"]
        assert len(delays) >= 3 and max(delays)<2, delays
        try:
            request("/api/sections",dict(os_entry=56,P=[1],Q=[1]))
        except urllib.error.HTTPError as error:
            assert error.code == 422
            assert "Neither is narrow" in json.load(error)["error"]
        else:
            raise AssertionError("Invalid pair was accepted")
        print("PASS: real OS56 result, %d concurrent health checks, max latency %.3fs; validation errors preserved" %
              (len(delays),max(delays)))
    finally:
        server.send_signal(2)
        try:
            server.wait(timeout=10)
        except subprocess.TimeoutExpired:
            server.kill()
            server.wait()
        if sys.exc_info()[0]:
            logs.seek(0)
            print(logs.read())
