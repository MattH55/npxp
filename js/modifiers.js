/* Modifiers: grouped list + detail panel with protocol parameters. */
renderNav("Modifiers");

const MOD_TYPE_LABEL = {
    dietary_metabolic: "Dietary / metabolic",
    thermal: "Thermal",
    hypoxic: "Hypoxic",
    mechanical_radiative: "Mechanical / radiative",
    other: "Other",
};
const MOD_TYPE_ICON = {
    dietary_metabolic: "fa-utensils",
    thermal: "fa-temperature-high",
    hypoxic: "fa-wind",
    mechanical_radiative: "fa-radiation",
    other: "fa-flask",
};

const state = { modifier_type: "", q: "" };

function protocolSummary(params) {
    const p = params || {};
    const bits = [];
    if (p.temperature_C) bits.push(`${p.temperature_C} °C`);
    if (p.concentration_mM) bits.push(`${p.concentration_mM} mM`);
    if (p.glucose_percent_of_standard) bits.push(`${p.glucose_percent_of_standard}% glucose`);
    if (p.o2_percent) bits.push(`${p.o2_percent}% O₂`);
    if (p.duration_hr) bits.push(`${p.duration_hr} h`);
    if (p.duration_min) bits.push(`${p.duration_min} min`);
    return bits.join(" · ") || "see detail";
}

async function loadList() {
    const status = document.getElementById("mod-status");
    status.textContent = "Loading…";
    try {
        const rows = await apiGet("/modifiers", state);
        status.textContent = `${rows.length} modifier(s)`;
        const groups = {};
        for (const m of rows) (groups[m.modifier_type] ||= []).push(m);
        document.getElementById("mod-groups").innerHTML = Object.entries(groups).map(([type, mods]) => `
          <section class="mb-6">
            <h2 class="text-white font-bold mb-2">
              <i class="fa-solid ${MOD_TYPE_ICON[type] || "fa-flask"} mr-2 text-emerald-400"></i>${MOD_TYPE_LABEL[type] || type}
              <span class="text-slate-500 text-sm font-normal ml-1">(${mods.length})</span>
            </h2>
            <div class="grid md:grid-cols-2 gap-3">
              ${mods.map((m) => `
                <div class="row-link bg-slate-800/60 border border-slate-700 rounded-xl p-4" data-id="${escapeHtml(m.modifier_id)}">
                  <div class="text-white font-semibold">${escapeHtml(m.agent)}</div>
                  <div class="text-slate-400 text-sm mt-1">${escapeHtml(protocolSummary(m.protocol_parameters))}</div>
                  <div class="flex items-center justify-between mt-2 text-xs text-slate-500">
                    <span>${m.signature_count} signature(s) · ${m.interaction_count} interaction(s)</span>
                    <span>${m.source_dataset_accession ? `<a class="cite-chip" target="_blank" rel="noopener" href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${escapeHtml(m.source_dataset_accession)}" onclick="event.stopPropagation()">GEO ${escapeHtml(m.source_dataset_accession)}</a>` : ""}</span>
                  </div>
                </div>`).join("")}
            </div>
          </section>`).join("") || '<div class="text-slate-500 text-sm">No modifiers match.</div>';
        document.querySelectorAll("#mod-groups .row-link").forEach((el) =>
            el.addEventListener("click", () => {
                location.hash = el.dataset.id;
                loadDetail(el.dataset.id);
            }));
    } catch (e) {
        status.textContent = `Failed to load modifiers: ${e.message}`;
    }
}

async function loadDetail(modifierId) {
    const box = document.getElementById("mod-detail");
    box.innerHTML = `<div class="text-slate-400 text-sm">Loading ${escapeHtml(modifierId)}…</div>`;
    try {
        const m = await apiGet(`/modifiers/${encodeURIComponent(modifierId)}`);
        const sigRows = (m.signature_scores || []).map((s) => `
          <tr class="border-b border-slate-700/60">
            <td class="py-2 pr-3 text-slate-300">${escapeHtml(s.cell_line_name || s.cell_line_id)}</td>
            <td class="py-2 pr-3 text-slate-300">${escapeHtml(s.signature_panel)}</td>
            <td class="py-2 pr-3 text-slate-300">${escapeHtml(s.directionality)}</td>
            <td class="py-2 pr-3 text-slate-300">${s.score === null ? '<span class="text-slate-500">n/a</span>' : escapeHtml(s.score)}</td>
            <td class="py-2 text-slate-400">${escapeHtml(s.notes)}</td>
          </tr>`).join("");
        const intRows = (m.interactions || []).map((it) => `
          <tr class="row-link border-b border-slate-700/60" onclick="location.href='interaction.html?id=${it.id}'">
            <td class="py-2 pr-3 text-white">${escapeHtml(it.drug_name)}${it.drug_lincs_signature_ref ? ' <i class="fa-solid fa-dna text-[0.6rem] text-emerald-400" title="LINCS-backed"></i>' : ""}</td>
            <td class="py-2 pr-3 text-slate-300">${escapeHtml(it.cell_line_name)}</td>
            <td class="py-2 pr-3">${typeBadge(it.interaction_type)}</td>
            <td class="py-2 pr-3">${tierBadge(it.evidence_tier)}</td>
            <td class="py-2">${citationChips(it.source_study)}</td>
          </tr>`).join("");
        box.innerHTML = `
          <div class="bg-slate-800/60 border border-slate-700 rounded-xl p-6">
            <h2 class="text-xl font-bold text-white mb-1">${escapeHtml(m.agent)}
              <span class="text-slate-500 text-sm font-normal ml-2">${escapeHtml(m.modifier_id)}</span></h2>
            <div class="text-slate-400 text-sm mb-4">
              ${escapeHtml(MOD_TYPE_LABEL[m.modifier_type] || m.modifier_type)} ·
              Source: ${citationChips([m.source_study])}
              ${m.source_dataset_accession ? citationChips(["geo:" + m.source_dataset_accession]) : ""}
            </div>
            <div class="text-xs text-slate-400 mb-4">
              Infrastructure requirement: <span class="cite-chip">${escapeHtml(m.infrastructure_requirement || "unassessed")}</span>
              ${m.estimated_relative_cost_note ? `<span class="ml-2">${escapeHtml(m.estimated_relative_cost_note)}</span>` : ""}
            </div>
            <h3 class="text-white font-semibold mb-2"><i class="fa-solid fa-sliders mr-2 text-emerald-400"></i>Protocol parameters</h3>
            <div class="grid sm:grid-cols-2 md:grid-cols-3 gap-3 mb-4">${protocolKv(m.protocol_parameters)}</div>
            <h3 class="text-white font-semibold mt-5"><i class="fa-solid fa-wave-square mr-2 text-emerald-400"></i>Stress signatures by cell line</h3>
            ${sigRows ? `<table class="w-full text-sm mt-2">
              <thead><tr class="text-left text-slate-400 border-b border-slate-700">
                <th class="py-2 pr-3">Cell line</th><th class="py-2 pr-3">Panel</th>
                <th class="py-2 pr-3">Directionality</th><th class="py-2 pr-3">Score</th><th class="py-2">Notes</th>
              </tr></thead><tbody>${sigRows}</tbody></table>` : '<div class="text-slate-500 text-sm">None recorded.</div>'}
            <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-circle-nodes mr-2 text-emerald-400"></i>Drug interactions</h3>
            ${intRows ? `<table class="w-full text-sm mt-2">
              <thead><tr class="text-left text-slate-400 border-b border-slate-700">
                <th class="py-2 pr-3">Drug</th><th class="py-2 pr-3">Cell line</th>
                <th class="py-2 pr-3">Type</th><th class="py-2 pr-3">Tier</th><th class="py-2">Citations</th>
              </tr></thead><tbody>${intRows}</tbody></table>` : '<div class="text-slate-500 text-sm">None recorded.</div>'}
          </div>`;
        box.scrollIntoView({ behavior: "smooth", block: "start" });
    } catch (e) {
        box.innerHTML = `<div class="text-rose-400 text-sm">${escapeHtml(e.message)}</div>`;
    }
}

let debounce;
document.getElementById("flt-type").addEventListener("change", () => {
    state.modifier_type = document.getElementById("flt-type").value;
    loadList();
});
document.getElementById("flt-q").addEventListener("input", () => {
    clearTimeout(debounce);
    debounce = setTimeout(() => {
        state.q = document.getElementById("flt-q").value.trim();
        loadList();
    }, 250);
});

loadList().then(() => {
    if (location.hash) loadDetail(decodeURIComponent(location.hash.slice(1)));
});
