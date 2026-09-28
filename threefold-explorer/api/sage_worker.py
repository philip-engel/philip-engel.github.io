"""Persistent Sage worker. All Sage calls execute on its main thread."""
from pathlib import Path
import json
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "engine"))
# Keep ordinary diagnostic prints out of the line-based response protocol.
PROTOCOL = sys.stdout
sys.stdout = sys.stderr
from explorer_api import API_VERSION, compute, describe_os_entry, log_transform_schema, section_and_linearization_schema
from local_model_database import LocalModelDatabase

DATABASE = LocalModelDatabase()

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
        return compute(body, verbose=False, database=DATABASE)
    raise KeyError("Unknown endpoint")


def send(value):
    print(json.dumps(value, separators=(",", ":")), file=PROTOCOL, flush=True)


if __name__ == "__main__":
    send(dict(status="ok", models=len(DATABASE.index["models"]),
              api_version=API_VERSION, scope="fiberwise-narrow"))
    for line in sys.stdin:
        started = time.monotonic()
        request = json.loads(line)
        path = request["path"]
        print("Sage request started: %s" % path, flush=True)
        try:
            response = dict(status=200, body=payload_for(path, request["body"]))
        except KeyError as error:
            response = dict(status=404 if str(error) == "'Unknown endpoint'" else 422,
                            body=dict(error=str(error).strip("'")))
        except (ValueError, TypeError, ArithmeticError, NotImplementedError) as error:
            response = dict(status=422, body=dict(error=str(error)))
        except Exception:
            traceback.print_exc()
            response = dict(status=500, body=dict(error="The Sage computation failed unexpectedly."))
        send(response)
        print("Sage request finished: %s, status %d, %.2fs" %
              (path, response["status"], time.monotonic()-started), flush=True)
