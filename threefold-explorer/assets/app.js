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
  let loadingPreset = false;
  let rationalSpherePresets;
  let lastRandomPresetIndex = -1;

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

  function stageError(id, error = null) {
    const box = $(id);
    box.textContent = error ? `Error: ${error.message}` : "";
    box.hidden = !error;
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
    if (body.api_version !== 2) {
      throw new Error("The Sage service is still running the previous version. Please try again after the update finishes.");
    }
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
      if (result.api_version !== 2) {
        setService("checking", "Sage service updating…");
        return false;
      }
      setService("online", `Sage online · ${result.models} local models`);
      return true;
    } catch {
      setService("offline", "Sage service unavailable");
      return false;
    }
  }

  function parseVector(value, name) {
    const pieces = value.trim().split(/[\s,]+/).filter(Boolean);
    if (!pieces.length && state.surface?.mordell_weil.tuple_length === 0) return [];
    if (!pieces.length) throw new Error(`${name} cannot be empty.`);
    return pieces.map((piece) => {
      if (!/^-?\d+$/.test(piece)) throw new Error(`${name} must contain integers separated by commas.`);
      return piece;
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
    const constraints = surface.fibers.filter((fiber) => fiber.Q_narrow_constraints.length)
      .map((fiber) => `<li><b>${fiber.type}:</b> ${fiber.narrowness_rule}</li>`);
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
    const narrowAt = pair.slots.filter((slot) => !slot.smooth).map((slot) =>
      `${slot.type}: ${slot.P_narrow && slot.Q_narrow ? "both" : slot.P_narrow ? "P" : "Q"}`);
    $("#pair-summary").textContent = `Pairing ⟨P,Q⟩ = ${pair.pairing}. Local narrowness satisfied (${narrowAt.join("; ")}).`;
    $("#pair-summary").classList.remove("pending", "error");
    $("#pair-summary").hidden = false;
    lock("#linearization-step", false);
    lock("#logs-step", true);
    $("#results").hidden = true;
  }

  function invalidatePairing() {
    stageError("#linearization-error");
    stageError("#logs-error");
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
    $("#pair-summary").textContent = `Error: ${error.message}`;
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
    stageError("#logs-error");
    state.logSchema = null;
    lock("#logs-step", true);
    $("#results").hidden = true;
    const currentWeights = weights();
    const total = currentWeights.reduce((sum, value) => sum + value, 0);
    $("#degree-total").textContent = String(total);
    const matches = state.pair && String(total) === state.pair.required_total_degree;
    $("#degree-total").style.color = matches ? "var(--green)" : "var(--ink)";
    let error = "";
    (state.pair?.slots || []).forEach((slot, index) => {
      const order = Math.abs(currentWeights[index] || 0);
      const match = /^I(\d+)$/.exec(slot.type);
      if (!order || !match) return;
      const components = Math.max(1, Number(match[1])) * order;
      if (components > 12) {
        error = "Error: The Mumford construction exceeds the allowed number of components (12).";
      } else if (order > 12 && !error) {
        error = "Error: The Mumford construction exceeds the allowed linearization order (12).";
      }
    });
    $("#linearization-error").textContent = error;
    $("#linearization-error").hidden = !error;
    $("#prepare-logs").disabled = Boolean(error);
  }

  function renderLogs(schema) {
    stageError("#logs-error");
    state.logSchema = schema;
    const ambient = $("#coordinates").value === "ambient";
    $("#log-row").innerHTML = schema.sites.map((site) => {
      const count = ambient ? 4 : site.coordinate_count;
      const orderHelp = site.arbitrary_denominators ?
        "Any torsion order is allowed at this smooth fiber." :
        `Allowed denominators divide ${site.reduction_order}.`;
      const mumfordHelp = site.type === "I0" && site.weight !== 0 ?
        " This slot has a Mumford filling. For a fractional twist, add a separate smooth slot with weight 0." : "";
      const basis = (ambient ? "Coordinates in (e₁,e₂,δ,c)." :
        `${site.coordinate_count}-dimensional invariant lattice.`) + ` ${orderHelp}${mumfordHelp}`;
      const columns = Array.from({length: site.coordinate_count}, (_, j) =>
        `v${j+1} = (${site.invariant_basis.map((row) => row[j]).join(", ")})`);
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
        ${ambient ? "" : `<details class="basis-help"><summary>Invariant basis in (e₁,e₂,δ,c)</summary>${columns.map((column) => `<div><code>${column}</code></div>`).join("")}</details>`}
      </article>`;
    }).join("");
    document.querySelectorAll(".log-mode").forEach((select) => select.addEventListener("change", () => {
      select.closest(".log-card").querySelector(".log-vector").hidden = select.value !== "vector";
    }));
    lock("#logs-step", false);
    $("#results").hidden = true;
  }

  async function prepareLogs() {
    notice("");
    stageError("#linearization-error");
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
      return `Resolved non-free good-reduction quotient${parameters.denominator ? ` by a cyclic group of order ${parameters.denominator}` : ""}${multiplicity ? `; fiber multiplicity ${multiplicity}` : ""}.`;
    }
    if (model.family === "star_orbit_quotient") {
      return `Minimal resolved quadratic quotient with Q permuting the components of the upstairs I${2 * Number(parameters.n || 0)} fiber${multiplicity ? `; fiber multiplicity ${multiplicity}` : ""}.`;
    }
    if (model.family === "star_semistable_quotient") {
      return `Minimal resolved quadratic quotient of the semistable I${2 * Number(parameters.n || 0)} model${multiplicity ? `; fiber multiplicity ${multiplicity}` : ""}.`;
    }
    if (model.family === "mumford") {
      const weight = Math.abs(Number(parameters.weight || 0));
      return `Semistable Mumford filling${weight ? ` with linearization order ${weight}` : ""} and A₂ tiling.`;
    }
    if (model.family === "original_plumbing") {
      if (model.component_orbits != null) {
        return `Original filling determined by O(P−O). Translation by Q has ${model.component_orbits} component orbit${model.component_orbits === 1 ? "" : "s"}, each of length ${model.component_orbit_length}; the reduced fiber is a wheel of ${model.component_orbits} component${model.component_orbits === 1 ? "" : "s"}.`;
      }
      return "Original filling determined by O(P−O), with no good-reduction substitution.";
    }
    if (model.family === "smooth_product") return multiplicity > 1 ?
      `Non-reduced fiber of multiplicity ${multiplicity}, with reduction a complex 2-torus.` :
      "Smooth T⁴ filling over a disk.";
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
        <span>${model.family === "smooth_product" && model.fiber_multiplicity > 1 ? "multiple torus" : model.family.replaceAll("_", " ")}</span>
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
    stageError("#logs-error");
    $("#results").hidden = true;
    try {
      const payload = currentBasePayload();
      payload.linearization_divisor = weights();
      payload.log_data = logData();
      payload.coordinates = $("#coordinates").value;
      const result = await post("/compute", payload);
      renderResult(result);
    } catch (error) { stageError("#logs-error", error); }
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
    cuspidal: {
      osEntry: 43, profile: "II", P: "1", Q: "2", weights: [0, 0, 1],
      vectors: { 0: "0, -1/4" },
      label: "The III* + II + I₁ example is ready. Press “Compute topology”.",
    },
    "os49-default": {
      osEntry: 49, profile: "default", P: "2", Q: "3", weights: [0, 0, 1, 0],
      vectors: { 1: "0, 1/3" },
      label: "OS49: I₂ + IV* + 2 I₁, with the full order twist at IV*. Press “Compute topology”.",
    },
    "os45-ii": {
      osEntry: 45, profile: "II", P: "8", Q: "1", weights: [0, 0, 1, 0],
      vectors: { 1: "0, 1/6" },
      label: "OS45: I₈ + II + 2 I₁, with the full order twist at II. Press “Compute topology”.",
    },
    "os47-ii": {
      osEntry: 47, profile: "II", P: "14", Q: "1", weights: [0, 0, 1, 0],
      vectors: { 3: "0, 1/6" },
      label: "OS47: I₇ + I₂ + I₁ + II, with the full order twist at II. Press “Compute topology”.",
    },
    "os47-iii": {
      osEntry: 47, profile: "III", P: "14", Q: "1", smoothSlots: 1,
      weights: [0, 1, 0, 0, 0], vectors: { 4: "0, 0, 0, 1" },
      label: "OS47: I₇ + III + 2 I₁, with one simple zero at I₁ and primitive integral clutching at the added smooth fiber. Press “Compute topology”.",
    },
    "os47-iii-p7": {
      osEntry: 47, profile: "III", P: "7", Q: "2", weights: [0, 1, 0, 0],
      vectors: { 3: "0, 1/4" },
      label: "OS47: P = 7G, Q = 2G, with the full order four twist at III. Press “Compute topology”.",
    },
    "os49-iii-p2": {
      osEntry: 49, profile: "III", P: "2", Q: "3", weights: [0, 0, 1],
      vectors: { 0: "0, 1/3" },
      label: "OS49: P = 2G, Q = 3G, with the full order three twist at IV* and the original filling at III. Press “Compute topology”.",
    },
    "os55-ii": {
      osEntry: 55, profile: "II", P: "20", Q: "1", weights: [0, 0, 0, 1],
      vectors: { 2: "0, 1/6" },
      label: "OS55: I₅ + I₄ + II + I₁, with the full order twist at II. Press “Compute topology”.",
    },
    "os56-iv": {
      osEntry: 56, profile: "IV", P: "10", Q: "3", weights: [0, 0, 0, 1],
      vectors: { 2: "0, 1/3" },
      label: "OS56: I₅ + I₂ + IV + I₁, with the full order twist at IV. Press “Compute topology”.",
    },
    "os56-iii-p15": {
      osEntry: 56, profile: "III", P: "15", Q: "2", weights: [0, 0, 0, 1],
      vectors: { 2: "0, 1/4" },
      label: "OS56: P = 15G, Q = 2G, with the full order four twist at III. Press “Compute topology”.",
    },
    "os56-iii-p30": {
      osEntry: 56, profile: "III", P: "30", Q: "1", smoothSlots: 1,
      weights: [0, 0, 0, 1, 0], vectors: { 4: "0, 0, 0, 1" },
      label: "OS56: P = 30G, Q = G, with original singular fillings and primitive integral clutching at the added smooth fiber. Press “Compute topology”.",
    },
    "os56-iv-p30": {
      osEntry: 56, profile: "IV", P: "30", Q: "1", smoothSlots: 1,
      weights: [0, 0, 0, 1, 0], vectors: { 4: "0, 0, 0, 1" },
      label: "OS56: P = 30G, Q = G, with the original filling at IV and primitive integral clutching at the added smooth fiber. Press “Compute topology”.",
    },
  };
  for (const [entry, multiple] of [[45, 8], [47, 14], [55, 20], [56, 30]]) {
    examplePresets[`os${entry}`] = {
      osEntry: entry, profile: "default", P: String(multiple), Q: "1",
      smoothSlots: 1, weights: [0, 0, 0, 0, 1, 0],
      vectors: {5: "0, 0, 0, 1"},
      label: `OS${entry}: Q = G, P = ${multiple}G, one simple zero at I₁, and primitive integral clutching at the added smooth fiber. Press “Compute topology”.`,
    };
  }

  async function loadExample(name) {
    await loadPreset(() => examplePresets[name]);
  }

  function randomPresetIndex(count, previous) {
    // Draw uniformly, omitting the previous choice after the first click.
    const omitPrevious = previous >= 0 && previous < count;
    const index = Math.floor(Math.random() * (count - (omitPrevious ? 1 : 0)));
    return omitPrevious && index >= previous ? index + 1 : index;
  }

  async function loadRandomExample() {
    let selectedIndex;
    const loaded = await loadPreset(async () => {
      if (!rationalSpherePresets) {
        const response = await fetch("assets/rational-sphere-presets.json");
        if (!response.ok) throw new Error("The random examples could not be loaded. Please try again.");
        const presets = await response.json();
        if (!Array.isArray(presets) || presets.length !== 99) {
          throw new Error("The random examples are unavailable. Please reload the page.");
        }
        rationalSpherePresets = presets;
      }
      selectedIndex = randomPresetIndex(rationalSpherePresets.length, lastRandomPresetIndex);
      const preset = rationalSpherePresets[selectedIndex];
      return {
        ...preset,
        label: `A random Q-homology sphere is ready (OS${preset.osEntry}, P = ${preset.P}, Q = ${preset.Q}). Press “Compute topology”.`,
      };
    });
    if (loaded) lastRandomPresetIndex = selectedIndex;
  }

  async function loadPreset(resolvePreset) {
    if (loadingPreset) return false;
    loadingPreset = true;
    clearTimeout(surfaceReloadTimer);
    const buttons = [$("#load-s6-example"), $("#load-rational-example")].filter(Boolean);
    buttons.forEach((button) => setBusy(button, true, "Loading…"));
    $("#compute").disabled = true;
    stageError("#linearization-error");
    stageError("#logs-error");
    try {
      const preset = await resolvePreset();
      $("#os-entry").value = String(preset.osEntry);
      await loadSurface(preset.profile);
      $("#p-vector").value = preset.P;
      $("#q-vector").value = preset.Q;
      state.smoothSlots = preset.smoothSlots || 0;
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
      return true;
    } catch (error) {
      notice(error.message, true);
      return false;
    } finally {
      loadingPreset = false;
      buttons.forEach((button) => setBusy(button, false));
      $("#compute").disabled = false;
    }
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
  $("#load-s6-example").addEventListener("click", () => loadExample($("#s6-example").value));
  $("#load-rational-example")?.addEventListener("click", loadRandomExample);
  $("#sections-form").addEventListener("submit", async (event) => {
    event.preventDefault(); notice("");
    try { await checkSections(); } catch (error) { showSectionError(error); }
  });
  $("#p-vector").addEventListener("input", invalidatePairing);
  $("#q-vector").addEventListener("input", invalidatePairing);
  $("#add-smooth").addEventListener("click", async () => {
    const savedWeights = weights();
    state.smoothSlots += 1;
    try {
      await checkSections();
      document.querySelectorAll(".weight-input").forEach((input, index) => {
        input.value = savedWeights[index] || 0;
      });
      updateDegree();
    } catch (error) { state.smoothSlots -= 1; stageError("#linearization-error", error); }
  });
  $("#prepare-logs").addEventListener("click", async () => {
    try { await prepareLogs(); } catch (error) { stageError("#linearization-error", error); }
  });
  $("#log-row").addEventListener("input", () => {
    stageError("#logs-error");
    $("#results").hidden = true;
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
