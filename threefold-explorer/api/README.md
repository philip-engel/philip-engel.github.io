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

Set `ALLOWED_ORIGINS` to a comma-separated list of exact HTTPS origins. The
default permits the GitHub Pages site and local development on port 8080.

```bash
docker build -t threefold-explorer-api .
docker run --rm -p 8000:8000 threefold-explorer-api
```
