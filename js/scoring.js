/* Synergy scoring tool (Prediction Layer, Step 1) — client for POST /api/scoring/synergy. */

renderNav("Scoring");

let AVAILABLE_MODELS = [];

function parseMatrix(text) {
    const lines = text.split(/\r?\n/).map((l) => l.trim()).filter((l) => l.length);
    if (lines.length < 2) throw new Error("Need a header row plus at least one modifier-level row.");
    const rows = lines.map((l) => l.split(/[,\t]/).map((c) => c.trim()));
    const width = rows[0].length;
    if (width < 3) throw new Error("Header row needs a label cell plus control (0) plus >=1 drug dose.");
    rows.forEach((r, i) => {
        if (r.length !== width) throw new Error(`Row ${i + 1} has ${r.length} cells; expected ${width}.`);
    });
    const num = (s, where) => {
        const v = Number(s);
        if (s === "" || Number.isNaN(v)) throw new Error(`Non-numeric value '${s}' at ${where}.`);
        return v;
    };
    const doses_b = rows[0].slice(1).map((s, j) => num(s, `header col ${j + 2}`));
    const doses_a = [];
    const viability = [];
    rows.slice(1).forEach((r, i) => {
        doses_a.push(num(r[0], `row ${i + 2}, first cell`));
        viability.push(r.slice(1).map((s, j) => num(s, `row ${i + 2}, col ${j + 2}`)));
    });
    return { doses_a, doses_b, viability };
}

function scoreColor(s) {
    if (s === null || s === undefined) return "background:#1e293b;color:#475569";
    const a = Math.min(1, Math.abs(s) / 0.5) * 0.85 + 0.1;
    return s >= 0
        ? `background:rgba(16,185,129,${a.toFixed(2)});color:#fff`
        : `background:rgba(56,189,248,${a.toFixed(2)});color:#fff`;
}

function renderResult(name, model, doses_a, doses_b) {
    const s = model.summary;
    const summary = s.n_cells === 0
        ? '<span class="text-slate-500 text-sm">No scorable combination cells.</span>'
        : `${typeBadge(s.classification)}
           <span class="text-slate-300 text-sm ml-2">mean ${s.mean.toFixed(3)} · max ${s.max.toFixed(3)} · min ${s.min.toFixed(3)} · ${s.n_cells} cells</span>`;
    const ciLine = model.meta && model.meta.mean_ci !== null && model.meta.mean_ci !== undefined
        ? `<span class="text-slate-400 text-sm ml-2">· mean CI ${model.meta.mean_ci.toFixed(3)}</span>` : "";
    const head = `<tr><th class="px-2 py-1 text-slate-500 text-xs font-mono">lvl\\dose</th>${doses_b.map((d) =>
        `<th class="px-2 py-1 text-slate-400 text-xs font-mono">${escapeHtml(d)}</th>`).join("")}</tr>`;
    const body = model.scores.map((row, i) =>
        `<tr><th class="px-2 py-1 text-slate-400 text-xs font-mono">${escapeHtml(doses_a[i])}</th>${row.map((v, j) => {
            const edge = i === 0 || j === 0;
            const label = v === null || v === undefined ? (edge ? "edge" : "n/a") : v.toFixed(3);
            return `<td class="px-2 py-1 text-xs font-mono text-center rounded" style="${edge ? "background:#0f172a;color:#334155" : scoreColor(v)}">${label}</td>`;
        }).join("")}</tr>`).join("");
    return `<div class="bg-slate-800/60 border border-slate-700 rounded-xl p-5 mb-5">
      <h2 class="text-white font-semibold mb-2">${escapeHtml(name.toUpperCase())}</h2>
      <div class="mb-3">${summary}${ciLine}</div>
      <table class="border-separate" style="border-spacing:3px">${head}${body}</table>
    </div>`;
}

async function runScoring() {
    const errEl = document.getElementById("input-error");
    errEl.textContent = "";
    let parsed;
    try {
        parsed = parseMatrix(document.getElementById("matrix-input").value);
    } catch (e) {
        errEl.textContent = e.message;
        return;
    }
    const models = AVAILABLE_MODELS.filter((m) => document.getElementById(`model-${m}`).checked);
    if (!models.length) { errEl.textContent = "Select at least one model."; return; }
    const scale = document.getElementById("scale").value;
    let data;
    if (typeof STATIC_DATA !== "undefined" && STATIC_DATA) {
        // Static build: same math runs in the browser (js/scoring-client.js).
        try {
            data = scoreMatrixClient(parsed.doses_a, parsed.doses_b, parsed.viability, models, scale);
        } catch (e) {
            errEl.textContent = e.message;
            return;
        }
    } else {
        const res = await fetch("/api/scoring/synergy", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ ...parsed, models, viability_scale: scale }),
        });
        if (!res.ok) {
            let detail = res.statusText;
            try { detail = (await res.json()).detail || detail; } catch (e) { /* non-JSON */ }
            errEl.textContent = `API ${res.status}: ${detail}`;
            return;
        }
        data = await res.json();
    }
    document.getElementById("results").innerHTML = Object.entries(data.models)
        .map(([name, model]) => renderResult(name, model, parsed.doses_a, parsed.doses_b)).join("");
}

async function init() {
    try { AVAILABLE_MODELS = (await apiGet("/scoring/models")).map((m) => m.name); }
    catch (e) { AVAILABLE_MODELS = ["bliss", "hsa", "loewe", "zip"]; }
    document.getElementById("model-checks").innerHTML = AVAILABLE_MODELS.map((m) =>
        `<label class="text-slate-300"><input type="checkbox" id="model-${m}" checked class="mr-1 accent-emerald-500">${m}</label>`).join("");
    try {
        const r = await apiGet("/prediction/readiness");
        document.getElementById("readiness").innerHTML = `
          <div class="border border-violet-500/40 bg-violet-500/10 rounded-lg p-3 text-sm text-violet-200">
            <i class="fa-solid fa-flask mr-1"></i>Step 2 (model-predicted interactions, Tier 2b) is
            <b>${r.ready ? "ready" : "disabled — insufficient ground truth"}</b>:
            ${r.tier_1_quantitative_count}/${r.required_quantitative} quantitative Tier 1 labels.
            Ship Step 1 labels first; predictions stay marked experimental until then.
          </div>`;
    } catch (e) { /* readiness endpoint optional */ }
    if (typeof STATIC_DATA !== "undefined" && STATIC_DATA) {
        document.getElementById("readiness").insertAdjacentHTML("beforeend",
            '<div class="text-slate-500 text-xs mt-2">Static build: scoring runs in your browser ' +
            "(same reference-model math as the server engine); matrices never leave this page.</div>");
    }
    document.getElementById("score-btn").addEventListener("click", runScoring);
}

init();
