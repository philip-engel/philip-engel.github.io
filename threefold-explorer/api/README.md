# Sage service

The service exposes the compact Threefold Explorer engine over JSON:

- `POST /api/os-entry`
- `POST /api/sections`
- `POST /api/log-schema`
- `POST /api/compute`
- `GET /health`

It uses the official `sagemath/sagemath:10.9` container and has no additional
Python dependencies. One service process handles requests serially because the
underlying Sage routines are not thread-safe. Run multiple container replicas
for concurrency.

The Git repository stores the database as `engine/local_model_database.tar.gz`
to avoid hundreds of tiny source files. The Docker build expands it into the
path expected by the engine. A local non-container run may either retain an
expanded development copy there or unpack the archive first.

The archive contains 865 reduced local cochain models and 450 shared exact
attachment formula tables. Ship the archive and engine modules together:
`symbolic_attachments.py` resolves the formula references and evaluates the
integral marking and clutching maps. A warm worker caches decoded records;
requests do not build local models or compile formula tables.

Weight-zero smooth I0 slots accept fractional log vectors of any torsion
order. They reuse the smooth torus record with an exact integral boundary
marking; the archive does not grow with the denominator. `/api/log-schema`
reports `arbitrary_denominators: true` and `allowable_denominators: null`
for these sites. Fractional twists on Mumford or other singular fibers remain
subject to the existing scope. Run `tests/smooth_logs_test.py` with Sage for
the local lattice identities and the single/coprime double-twist regressions.

Run `sage -python tests/worker_bundle_test.py` to check the shipped archive in
an isolated temporary directory through the persistent worker protocol. Run
`sage -python tests/smoke_test.py` to check the expanded database, and
`sage -python tests/service_test.py` for HTTP integration where local listening
sockets are permitted.

Set `ALLOWED_ORIGINS` to a comma-separated list of exact HTTPS origins. The
default permits the GitHub Pages site and local development on port 8080.

```bash
docker build -t threefold-explorer-api .
docker run --rm -p 8000:8000 threefold-explorer-api
```

The OS collision lookup is `engine/collision_profiles.json`. Keep it beside
`os_monodromy.py`; the Docker build and isolated worker test include it. It
contains 215 alternate profiles, covering 279 fiber configurations together
with the 74 defaults. Run `tests/collision_profiles_test.py` with Sage to check
all 289 marked profiles.
