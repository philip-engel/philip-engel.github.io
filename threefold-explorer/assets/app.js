(() => {
  "use strict";

  const $ = (selector) => document.querySelector(selector);
  const state = {
    apiBase: localStorage.getItem("threefoldExplorerApi") ||
      window.THREEFOLD_EXPLORER_CONFIG?.apiBase || "",
    surface: null,
    pair: null,
    logSchema: null,
    smoothSlots: 0,
  };
  let surfaceReloadTimer;

  function setService(kind, label) {
    $("#service-state").dataset.state = kind;
    $("#service-label").textContent = label;
  }

  function notice(message, error = false) {
    const box = $("#notice");
    box.textContent = message;
    box.classList.toggle("error", error);
    box.hidden = !message;
  }

  function lock(id, locked) {
    const element = $(id);
    element.classList.toggle("locked", locked);
    element.classList.toggle("active", !locked);
    element.setAttribute("aria-disabled", String(locked));
  }

  function apiUrl(path) {
    return `${state.apiBase.replace(/\/$/, "")}${path}`;
  }

  async function post(path, payload) {
    if (!state.apiBase) throw new Error("Enter the Sage service address under “Sage service”.");
    const response = await fetch(apiUrl(path), {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    let body;
    try { body = await response.json(); }
    catch { throw new Error("The Sage service returned an unreadable response."); }
    if (!response.ok) throw new Error(body.error || `Request failed (${response.status}).`);
    return body;
  }

  async function checkService() {
    if (!state.apiBase) {
      setService("offline", "Sage service not configured");
      return false;
    }
    try {
      const root = state.apiBase.replace(/\/api\/?$/, "");
      const response = await fetch(`${root}/health`, { cache: "no-store" });
      const result = await response.json();
      if (!response.ok || result.status !== "ok") throw new Error();
      setService("online", `Sage online · ${result.models} local models`);
      return true;
    } catch {
      setService("offline", "Sage service unavailable");
      return false;
    }
  }

  function parseVector(value, name) {
    const pieces = value.trim().split(/[\s,]+/).filter(Boolean);
    if (!pieces.length) throw new Error(`${name} cannot be empty.`);
    return pieces.map((piece) => {
      if (!/^-?\d+$/.test(piece)) throw new Error(`${name} must contain integers separated by commas.`);
      return Number(piece);
    });
  }

  function parseRationals(value, expected, name) {
    const pieces = value.trim().split(/[\s,]+/).filter(Boolean);
    if (pieces.length !== expected) throw new Error(`${name} requires ${expected} entries.`);
    for (const piece of pieces) {
      if (!/^-?\d+(?:\/-?\d+)?$/.test(piece) || /\/0$/.test(piece)) {
        throw new Error(`${name} must use exact integers or fractions such as 1/2.`);
      }
    }
    return pieces;
  }

  function currentBasePayload() {
    if (!state.surface) throw new Error("Load an OS entry first.");
    return {
      os_entry: state.surface.os_entry,
      profile: state.surface.profile,
      P: parseVector($("#p-vector").value, "P"),
      Q: parseVector($("#q-vector").value, "Q"),
    };
  }

  function setBusy(button, busy, busyText) {
    if (!button.dataset.label) button.dataset.label = button.textContent;
    button.disabled = busy;
    button.textContent = busy ? busyText : button.dataset.label;
  }

  function renderHeightPairing(surface) {
    const rank = surface.mordell_weil.rank;
    const matrix = surface.mordell_weil.height_matrix;
    const panel = $("#height-panel");
    if (!rank) {
      $("#height-matrix").innerHTML = "";
      $("#height-formula").textContent = "The height pairing vanishes because MW has no free part.";
      panel.hidden = false;
      return;
    }
    $("#height-matrix").innerHTML = `<table class="height-matrix" aria-label="Height pairing matrix"><tbody>${matrix.map((row) =>
      `<tr>${row.map((entry) => `<td>${entry}</td>`).join("")}</tr>`
    ).join("")}</tbody></table>`;
    $("#height-formula").textContent = "⟨P,Q⟩ = pᵀHq on the free coordinates; torsion coordinates do not contribute.";
    panel.hidden = false;
  }

  function renderSurface(surface) {
    state.surface = surface;
    state.pair = null;
    state.logSchema = null;
    state.smoothSlots = 0;
    $("#profile").innerHTML = surface.collision_profiles.map((profile) =>
      `<option value="${profile.profile}" ${profile.selected ? "selected" : ""}>${profile.label}</option>`
    ).join("");
    $("#fiber-ribbon").innerHTML = surface.fibers.map((fiber) =>
      `<div class="fiber-chip"><small>fiber ${fiber.index}</small>${fiber.type}</div>`
    ).join("");
    $("#mw-label").textContent = `MW = ${surface.mordell_weil.label}`;
    $("#tuple-label").textContent = `${surface.mordell_weil.tuple_length} section coordinate${surface.mordell_weil.tuple_length === 1 ? "" : "s"}`;
    renderHeightPairing(surface);
    $("#surface-summary").hidden = false;
    const zeros = Array(surface.mordell_weil.tuple_length).fill(0).join(", ");
    $("#p-vector").placeholder = zeros;
    $("#q-vector").placeholder = zeros;
    const constraints = surface.fibers.flatMap((fiber) =>
      fiber.Q_narrow_constraints.map((item) => `<li><b>${fiber.type}:</b> ${item.display}</li>`)
    );
    $("#narrow-constraints").innerHTML = constraints.length
      ? `<ul class="constraint-list">${constraints.join("")}</ul>`
      : "<p>Every section is locally narrow at the displayed fibers.</p>";
    lock("#sections-step", false);
    lock("#linearization-step", true);
    lock("#logs-step", true);
    $("#results").hidden = true;
  }

  async function loadSurface(profileOverride) {
    const button = $("#surface-form .primary");
    setBusy(button, true, "Loading…");
    notice("");
    try {
      const osEntry = Number($("#os-entry").value);
      const entryChanged = state.surface && state.surface.os_entry !== osEntry;
      const surface = await post("/os-entry", {
        os_entry: osEntry,
        profile: profileOverride ?? (entryChanged ? "default" : $("#profile").value),
      });
      renderSurface(surface);
      return surface;
    } catch (error) {
      notice(error.message, true);
      throw error;
    } finally { setBusy(button, false); }
  }

  function renderLinearization(pair) {
    state.pair = pair;
    $("#degree-required").textContent = pair.required_total_degree;
    $("#linearization-row").innerHTML = pair.slots.map((slot) => `
      <label class="slot-card">
        <strong>${slot.type}</strong>
        <small>${slot.smooth ? "added smooth fiber" : `fiber ${slot.index}`} · ${slot.explanation}</small>
        <input class="weight-input" type="number" step="1" value="0" data-index="${slot.index - 1}" ${slot.linearization_allowed ? "" : "disabled"}>
      </label>`).join("");
    document.querySelectorAll(".weight-input").forEach((input) => input.addEventListener("input", updateDegree));
    updateDegree();
    $("#pair-summary").innerHTML = `Pairing ⟨P,Q⟩ = <b>${pair.pairing}</b>. Q satisfies every local narrowness condition.`;
    $("#pair-summary").classList.remove("pending", "error");
    $("#pair-summary").hidden = false;
    lock("#linearization-step", false);
    lock("#logs-step", true);
    $("#results").hidden = true;
  }

  function invalidatePairing() {
    if (!state.surface) return;
    state.pair = null;
    state.logSchema = null;
    $("#pair-summary").textContent = "P or Q changed. Press “Check sections” to recompute ⟨P,Q⟩.";
    $("#pair-summary").classList.remove("error");
    $("#pair-summary").classList.add("pending");
    $("#pair-summary").hidden = false;
    $("#degree-required").textContent = "—";
    $("#linearization-row").innerHTML = "";
    $("#log-row").innerHTML = "";
    lock("#linearization-step", true);
    lock("#logs-step", true);
    $("#results").hidden = true;
  }

  function showSectionError(error) {
    state.pair = null;
    state.logSchema = null;
    const narrownessFailure = /narrow/i.test(error.message);
    $("#pair-summary").textContent = narrownessFailure
      ? "Error: The current topology computation requires Q to be globally narrow."
      : error.message;
    $("#pair-summary").classList.remove("pending");
    $("#pair-summary").classList.add("error");
    $("#pair-summary").hidden = false;
    $("#degree-required").textContent = "—";
    $("#linearization-row").innerHTML = "";
    $("#log-row").innerHTML = "";
    lock("#linearization-step", true);
    lock("#logs-step", true);
    $("#results").hidden = true;
  }

  async function checkSections() {
    const payload = currentBasePayload();
    payload.smooth_slots = state.smoothSlots;
    const pair = await post("/sections", payload);
    renderLinearization(pair);
    return pair;
  }

  function weights() {
    return [...document.querySelectorAll(".weight-input")].map((input) => Number(input.value || 0));
  }

  function updateDegree() {
    const total = weights().reduce((sum, value) => sum + value, 0);
    $("#degree-total").textContent = String(total);
    const matches = state.pair && String(total) === state.pair.required_total_degree;
    $("#degree-total").style.color = matches ? "var(--green)" : "var(--ink)";
  }

  function renderLogs(schema) {
    state.logSchema = schema;
    const ambient = $("#coordinates").value === "ambient";
    $("#log-row").innerHTML = schema.sites.map((site) => {
      const count = ambient ? 4 : site.coordinate_count;
      const basis = ambient ? "Coordinates in (e₁,e₂,δ,c)." :
        `${site.coordinate_count}-dimensional invariant lattice; allowed denominators divide ${site.reduction_order}.`;
      return `<article class="log-card" data-index="${site.index - 1}" data-count="${count}">
        <h3>${site.type} <small>· fiber ${site.index}</small></h3>
        <label>Filling choice
          <select class="log-mode">
            <option value="none">None · ${site.none_meaning}</option>
            <option value="zero">Explicit zero · ${site.zero_meaning}</option>
            <option value="vector">Enter a log vector</option>
          </select>
        </label>
        <label class="log-vector" hidden>Exact vector (${count} entries)
          <input type="text" autocomplete="off" placeholder="${Array(count).fill("0").join(", ")}">
        </label>
        <p class="log-help">${basis}</p>
      </article>`;
    }).join("");
    document.querySelectorAll(".log-mode").forEach((select) => select.addEventListener("change", () => {
      select.closest(".log-card").querySelector(".log-vector").hidden = select.value !== "vector";
    }));
    lock("#logs-step", false);
    $("#results").hidden = true;
  }

  async function prepareLogs() {
    const payload = currentBasePayload();
    payload.linearization_divisor = weights();
    const schema = await post("/log-schema", payload);
    renderLogs(schema);
    $("#logs-step").scrollIntoView({ behavior: "smooth", block: "start" });
    return schema;
  }

  function logData() {
    return [...document.querySelectorAll(".log-card")].map((card, index) => {
      const mode = card.querySelector(".log-mode").value;
      const count = Number(card.dataset.count);
      if (mode === "none") return null;
      if (mode === "zero") return Array(count).fill(0);
      return parseRationals(card.querySelector("input").value, count, `Fiber ${index + 1} log vector`);
    });
  }

  function localGeometryDescription(model) {
    const parameters = model.parameters || {};
    const multiplicity = model.fiber_multiplicity ?? parameters.denominator;
    if (model.family === "finite_quotient" && parameters.resolution === "free") {
      return `Multiple fiber${multiplicity ? ` of multiplicity ${multiplicity}` : ""}; its reduction is a bielliptic surface obtained as a free cyclic quotient of the good-reduction abelian surface.`;
    }
    if (model.family === "finite_quotient") {
      return `Resolved non-free good-reduction quotient${multiplicity ? ` by a cyclic group of order ${multiplicity}` : ""}.`;
    }
    if (model.family === "star_semistable_quotient") {
      return `Minimal resolved quadratic quotient of the semistable I${2 * Number(parameters.n || 0)} model${multiplicity ? `; fiber multiplicity ${multiplicity}` : ""}.`;
    }
    if (model.family === "mumford") {
      const weight = Math.abs(Number(parameters.weight || 0));
      return `Semistable Mumford filling${weight ? ` with linearization order ${weight}` : ""} and A₂ tiling.`;
    }
    if (model.family === "original_plumbing") {
      return "Original filling determined by O(P−O), with no good-reduction substitution.";
    }
    if (model.family === "smooth_product") return "Smooth T⁴ filling over a disk.";
    return model.geometry || "";
  }

  function launchConfetti() {
    if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) return;
    document.querySelector(".confetti-layer")?.remove();
    const colors = ["#2854bd", "#3976e1", "#28785c", "#f1b83b", "#d85d6f", "#8a62cc"];
    const layer = document.createElement("div");
    layer.className = "confetti-layer";
    layer.setAttribute("aria-hidden", "true");
    for (let index = 0; index < 64; index += 1) {
      const piece = document.createElement("i");
      piece.style.setProperty("--x", `${Math.random() * 100}vw`);
      piece.style.setProperty("--drift", `${(Math.random() - .5) * 28}vw`);
      piece.style.setProperty("--spin", `${540 + Math.random() * 720}deg`);
      piece.style.setProperty("--delay", `${Math.random() * .45}s`);
      piece.style.setProperty("--duration", `${1.8 + Math.random() * .8}s`);
      piece.style.background = colors[index % colors.length];
      layer.appendChild(piece);
    }
    document.body.appendChild(layer);
    window.setTimeout(() => layer.remove(), 3200);
  }

  function renderResult(result) {
    $("#result-title").textContent = result.integral_homology_sphere ? "Integral homology sphere" : "Computed threefold";
    $("#sphere-badge").hidden = !result.S6_for_supplied_smooth_model;
    $("#pi1-result").textContent = result.fundamental_group.description;
    $("#h1-result").textContent = `Abelianization: ${result.fundamental_group.abelianization}`;
    $("#euler-result").textContent = String(result.euler_characteristic);
    $("#cohomology").innerHTML = result.cohomology.map((group) =>
      `<div class="cohomology-group"><small>H<sup>${group.degree}</sup></small><strong>${group.label}</strong></div>`
    ).join("");
    $("#local-models").innerHTML = result.local_models.map((model) =>
      `<div class="model-row"><span>${model.index}</span><b>${model.type}</b><div class="model-copy">
        <span>${model.family.replaceAll("_", " ")}</span>
        <small>${localGeometryDescription(model)}</small>
      </div></div>`
    ).join("");
    $("#results").hidden = false;
    $("#results").scrollIntoView({ behavior: "smooth", block: "start" });
    if (result.S6_for_supplied_smooth_model) launchConfetti();
  }

  async function runComputation() {
    const button = $("#compute");
    setBusy(button, true, "Computing with Sage…");
    notice("");
    try {
      const payload = currentBasePayload();
      payload.linearization_divisor = weights();
      payload.log_data = logData();
      payload.coordinates = $("#coordinates").value;
      const result = await post("/compute", payload);
      renderResult(result);
    } catch (error) { notice(error.message, true); }
    finally { setBusy(button, false); }
  }

  const examplePresets = {
    original: {
      osEntry: 49, profile: "III", P: "1", Q: "6", weights: [0, 0, 1],
      vectors: {
        0: "0, 1/3",
        1: "0, 1/4",
      },
      label: "The IV* + III + I₁ example is ready. Press “Compute topology”.",
    },
    split: {
      osEntry: 43, profile: "default", P: "1", Q: "2", weights: [0, 0, 0, 1],
      vectors: { 0: "0, -1/4" },
      label: "The III* + I₁ + I₁ + I₁ example is ready. Press “Compute topology”.",
    },
  };

  async function loadExample(name) {
    const preset = examplePresets[name];
    try {
      $("#os-entry").value = String(preset.osEntry);
      await loadSurface(preset.profile);
      $("#p-vector").value = preset.P;
      $("#q-vector").value = preset.Q;
      await checkSections();
      const fields = document.querySelectorAll(".weight-input");
      preset.weights.forEach((value, index) => { fields[index].value = value; });
      updateDegree();
      $("#coordinates").value = "invariant";
      await prepareLogs();
      Object.entries(preset.vectors).forEach(([index, vector]) => {
        const card = document.querySelector(`.log-card[data-index="${index}"]`);
        const mode = card.querySelector(".log-mode");
        mode.value = "vector";
        mode.dispatchEvent(new Event("change"));
        card.querySelector("input").value = vector;
      });
      notice(preset.label);
    } catch (error) { notice(error.message, true); }
  }

  $("#surface-form").addEventListener("submit", async (event) => {
    event.preventDefault();
    try { await loadSurface(); } catch { /* displayed above */ }
  });
  $("#os-entry").addEventListener("input", () => {
    clearTimeout(surfaceReloadTimer);
    surfaceReloadTimer = setTimeout(() => loadSurface("default").catch(() => {}), 300);
  });
  $("#profile").addEventListener("change", () => loadSurface().catch(() => {}));
  $("#load-original-example").addEventListener("click", () => loadExample("original"));
  $("#load-split-example").addEventListener("click", () => loadExample("split"));
  $("#sections-form").addEventListener("submit", async (event) => {
    event.preventDefault(); notice("");
    try { await checkSections(); } catch (error) { showSectionError(error); }
  });
  $("#p-vector").addEventListener("input", invalidatePairing);
  $("#q-vector").addEventListener("input", invalidatePairing);
  $("#add-smooth").addEventListener("click", async () => {
    state.smoothSlots += 1;
    try { await checkSections(); } catch (error) { state.smoothSlots -= 1; notice(error.message, true); }
  });
  $("#prepare-logs").addEventListener("click", async () => {
    try { await prepareLogs(); } catch (error) { notice(error.message, true); }
  });
  $("#coordinates").addEventListener("change", () => state.logSchema && renderLogs(state.logSchema));
  $("#compute").addEventListener("click", runComputation);
  $("#api-settings").addEventListener("click", () => {
    const panel = $("#settings-panel");
    panel.hidden = !panel.hidden;
    $("#api-settings").setAttribute("aria-expanded", String(!panel.hidden));
    $("#api-url").value = state.apiBase;
  });
  $("#save-api").addEventListener("click", async () => {
    state.apiBase = $("#api-url").value.trim().replace(/\/$/, "");
    localStorage.setItem("threefoldExplorerApi", state.apiBase);
    const online = await checkService();
    notice(online ? "Sage service connected." : "The Sage service could not be reached.", !online);
  });

  checkService();
})();
