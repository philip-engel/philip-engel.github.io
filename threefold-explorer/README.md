# Threefold Explorer web application

This directory contains the public interface and the containerized Sage service.

## Structure

- `index.html` and `assets/`: static GitHub Pages frontend;
- `api/`: SageMath 10.9 service, compressed 865-model read-only database with complete bounded Mumford and I_n* lookup tables under the fiberwise P-or-Q narrowness condition;
- `docker-compose.yml`: local full-stack launch;
- `docs/`: public scope and runtime documentation.

GitHub Pages serves the frontend at `/threefold-explorer/`. The frontend reads
the Sage service address from `assets/config.js`; visitors may also supply an
address through the service panel, which stores it locally in their browser.

## Local development

Start the Sage service directly on this Mac:

```bash
cd threefold-explorer/api
PORT=8000 sage -python service.py
```

In another terminal, serve the static site from the repository root:

```bash
python3 -m http.server 8080
```

Then open `http://127.0.0.1:8080/threefold-explorer/`.

For a container runtime, `docker compose up --build` starts the API. The static
frontend may still be served by any ordinary file server.

## Deployment

The GitHub Pages repository hosts the static frontend. `render.yaml` describes a
separate Docker web service built from `api/Dockerfile`. After that service is
created, put its HTTPS `/api` address in `assets/config.js` and keep
`ALLOWED_ORIGINS=https://philip-engel.github.io` on the service.

The HTTP process stays responsive while one persistent Sage subprocess handles computations serially. Health checks do not wait for a Sage calculation. Run the transport integration check with Sage’s Python and `api/tests/service_test.py`.

## Licensing

The program code is available under the [MIT License](LICENSE-CODE.txt). The
explanatory documentation is available under
[Creative Commons Attribution 4.0 International](LICENSE-DOCUMENTATION.md).
The tabulated local-model database and unpublished mathematical derivation
material are excluded from those licenses; see [COPYRIGHT.md](COPYRIGHT.md) for
the precise scope.
