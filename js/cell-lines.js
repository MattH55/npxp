/* Cell lines: filterable table + detail panel. */
renderNav("Cell Lines");

const state = { q: "", tissue: "", mutation: "" };

async function loadList() {
    const status = document.getElementById("cl-status");
    status.textContent = "Loading…";
    try {
        const rows = await apiGet("/cell_lines", state);
        status.textContent = `${rows.length} cell line(s)`;
        const tbody = document.getElementById("cl-tbody");
        tbody.innerHTML = rows.map((cl) => `
          <tr class="row-link border-b border-slate-700/60" data-id="${escapeHtml(cl.cell_line_id)}">
            <td class="px-4 py-3 text-white font-semibold">${escapeHtml(cl.name)}</td>
            <td class="px-4 py-3 text-slate-300">${escapeHtml(cl.tissue_origin)}</td>
            <td class="px-4 py-3 text-slate-400">${escapeHtml(cl.cancer_subtype)}</td>
            <td class="px-4 py-3">${mutationChips(cl.key_mutations)}</td>
            <td class="px-4 py-3 text-center text-slate-300">${cl.signature_count}</td>
            <td class="px-4 py-3 text-center text-slate-300">${cl.interaction_count}</td>
          </tr>`).join("") ||
          `<tr><td colspan="6" class="px-4 py-6 text-center text-slate-500">No cell lines match.</td></tr>`;
        tbody.querySelectorAll("tr.row-link").forEach((tr) =>
            tr.addEventListener("click", () => {
                location.hash = tr.dataset.id;
                loadDetail(tr.dataset.id);
            }));
    } catch (e) {
        status.textContent = `Failed to load cell lines: ${e.message}`;
    }
}

function interactionTable(interactions) {
    if (!interactions.length) return '<div class="text-slate-500 text-sm">No curated interactions.</div>';
    return `<table class="w-full text-sm mt-2">
      <thead><tr class="text-left text-slate-400 border-b border-slate-700">
        <th class="py-2 pr-3">Drug</th><th class="py-2 pr-3">Modifier</th>
        <th class="py-2 pr-3">Type</th><th class="py-2 pr-3">Tier</th><th class="py-2">Citations</th>
      </tr></thead>
      <tbody>${interactions.map((it) => `
        <tr class="row-link border-b border-slate-700/60" onclick="location.href='interaction.html?id=${it.id}'">
          <td class="py-2 pr-3 text-white">${escapeHtml(it.drug_name)}</td>
          <td class="py-2 pr-3 text-slate-300">${escapeHtml(it.modifier_agent)}</td>
          <td class="py-2 pr-3">${typeBadge(it.interaction_type)}</td>
          <td class="py-2 pr-3">${tierBadge(it.evidence_tier)}</td>
          <td class="py-2">${citationChips(it.source_study)}</td>
        </tr>`).join("")}
      </tbody></table>`;
}

function signatureTable(scores) {
    if (!scores.length) return '<div class="text-slate-500 text-sm">No stress signatures recorded.</div>';
    return `<table class="w-full text-sm mt-2">
      <thead><tr class="text-left text-slate-400 border-b border-slate-700">
        <th class="py-2 pr-3">Modifier</th><th class="py-2 pr-3">Panel</th>
        <th class="py-2 pr-3">Directionality</th><th class="py-2 pr-3">Score</th><th class="py-2">Notes</th>
      </tr></thead>
      <tbody>${scores.map((s) => `
        <tr class="border-b border-slate-700/60">
          <td class="py-2 pr-3 text-slate-300">${escapeHtml(s.modifier_agent || s.modifier_id)}</td>
          <td class="py-2 pr-3 text-slate-300">${escapeHtml(s.signature_panel)}</td>
          <td class="py-2 pr-3 text-slate-300">${escapeHtml(s.directionality)}</td>
          <td class="py-2 pr-3 text-slate-300">${s.score === null ? '<span class="text-slate-500">n/a*</span>' : escapeHtml(s.score)}</td>
          <td class="py-2 text-slate-400">${escapeHtml(s.notes)}</td>
        </tr>`).join("")}
      </tbody></table>
      <div class="text-slate-500 text-xs mt-1">* Quantitative ssGSEA/GSVA scores are backfilled by the ingestion pipeline; directionality is as reported in the cited source.</div>`;
}

async function loadDetail(cellLineId) {
    const box = document.getElementById("cl-detail");
    box.innerHTML = `<div class="text-slate-400 text-sm">Loading ${escapeHtml(cellLineId)}…</div>`;
    try {
        const cl = await apiGet(`/cell_lines/${encodeURIComponent(cellLineId)}`);
        box.innerHTML = `
          <div class="bg-slate-800/60 border border-slate-700 rounded-xl p-6">
            <h2 class="text-xl font-bold text-white mb-1">${escapeHtml(cl.name)}
              <span class="text-slate-500 text-sm font-normal ml-2">${escapeHtml(cl.cell_line_id)}</span></h2>
            <div class="text-slate-400 text-sm mb-4">${escapeHtml(cl.tissue_origin)} · ${escapeHtml(cl.cancer_subtype)} · Source: ${escapeHtml(cl.source)}</div>
            <div class="mb-4"><span class="kv-key mr-2">Key mutations</span>${mutationChips(cl.key_mutations)}</div>
            <h3 class="text-white font-semibold mt-5"><i class="fa-solid fa-wave-square mr-2 text-emerald-400"></i>Stress signatures</h3>
            ${signatureTable(cl.signature_scores)}
            <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-circle-nodes mr-2 text-emerald-400"></i>Drug × modifier interactions</h3>
            ${interactionTable(cl.interactions)}
          </div>`;
        box.scrollIntoView({ behavior: "smooth", block: "start" });
    } catch (e) {
        box.innerHTML = `<div class="text-rose-400 text-sm">${escapeHtml(e.message)}</div>`;
    }
}

let debounce;
["flt-q", "flt-tissue", "flt-mutation"].forEach((id) => {
    document.getElementById(id).addEventListener("input", () => {
        clearTimeout(debounce);
        debounce = setTimeout(() => {
            state.q = document.getElementById("flt-q").value.trim();
            state.tissue = document.getElementById("flt-tissue").value.trim();
            state.mutation = document.getElementById("flt-mutation").value.trim();
            loadList();
        }, 250);
    });
});
document.getElementById("flt-clear").addEventListener("click", () => {
    ["flt-q", "flt-tissue", "flt-mutation"].forEach((id) => (document.getElementById(id).value = ""));
    state.q = state.tissue = state.mutation = "";
    loadList();
});

loadList().then(() => {
    if (location.hash) loadDetail(decodeURIComponent(location.hash.slice(1)));
});
