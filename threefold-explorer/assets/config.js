window.THREEFOLD_EXPLORER_CONFIG = {
  apiBase: ["localhost", "127.0.0.1"].includes(window.location.hostname)
    ? "http://127.0.0.1:8000/api"
    : "https://threefold-explorer-api.onrender.com/api"
};
