/* Mechanistic Bridging candidates (build spec v5) + Prioritization Schema
   (build spec v6). The prioritized-worklist page: sorts by priority_score
   by default, with a published methodology panel fetched from the same
   endpoint the score is computed by (never a hardcoded second copy of the
   weights). Filtering/sorting is delegated to apiGet's params so the same
   code path works against the live API and the static JSON snapshot. */
renderNav("Candidates");

const state = { drug_id: "", tissue_origin: "", low_resource_relevance: "", sort: "priority_score" };

function candidateCard(it) {
    const rm = it.resistance_mechanism || {};
    return `
      <div class="row-link bg-slate-800/60 border border-violet-500/30 rounded-xl p-4"
           onclick="location.href='interaction.html?id=${it.id}'">
        <div class="flex flex-wrap items-center justify-between gap-2 mb-2">
          <div class="flex flex-wrap items-center gap-2">
            ${typeBadge(it.interaction_type)} ${tierBadge(it.evidence_tier)}
            <span class="cite-chip">model-generated</span>
            ${lowResourceBadge(it.low_resource_relevance)}
          </div>
          <div class="text-right">
            <div class="text-white font-bold text-lg leading-none">${it.priority_score}</div>
            <div class="text-slate-500 text-[0.65rem]">priority score</div>
          </div>
        </div>
        <div class="text-white font-semibold">
          ${escapeHtml(it.drug_name)}${it.drug_lincs_signature_ref ? ' <i class="fa-solid fa-dna text-[0.6rem] text-emerald-400" title="LINCS-backed"></i>' : ""}
          <span class="text-slate-500 font-normal">+</span> ${escapeHtml(it.modifier_agent)}
          <span class="text-slate-500 font-normal">in</span> ${escapeHtml(it.cell_line_name)}
        </div>
        <div class="text-slate-400 text-xs mt-1">${escapeHtml(it.tissue_origin || "")}</div>
        <div class="text-slate-300 text-sm mt-2">
          <span class="kv-key">Resistance pathway matched:</span> ${escapeHtml(rm.pathway_or_gene || "")}
        </div>
        <div class="text-slate-400 text-sm mt-1">${escapeHtml(it.curator_notes)}</div>
        <div class="mt-2">${citationChips(it.source_study)}</div>
      </div>`;
}

async function render() {
    const status = document.getElementById("cand-status");
    status.textContent = "Loading…";
    try {
        const params = { drug_id: state.drug_id, tissue_origin: state.tissue_origin, sort: state.sort };
        if (state.low_resource_relevance) params.low_resource_relevance = true;
        const rows = await apiGet("/interactions/candidates", params);
        status.textContent = `${rows.length} candidate(s)`;
        document.getElementById("cand-list").innerHTML =
            rows.map(candidateCard).join("") ||
            '<div class="text-slate-500 text-sm">No candidates match.</div>';
    } catch (e) {
        status.textContent = `Failed to load candidates: ${e.message}`;
    }
}

async function loadMethodology() {
    const box = document.getElementById("methodology");
    try {
        const meta = await apiGet("/prioritization/methodology");
        box.innerHTML = `
          <p class="mb-2">${escapeHtml(meta.formula)}</p>
          <table class="w-full text-xs mt-2">
            <thead><tr class="text-left text-slate-500 border-b border-slate-700">
              <th class="py-1 pr-3">Component</th><th class="py-1 pr-3">Weight</th><th class="py-1">What it measures</th>
            </tr></thead>
            <tbody>${Object.entries(meta.weights).map(([k, w]) => `
              <tr class="border-b border-slate-700/60">
                <td class="py-1 pr-3 text-slate-200">${escapeHtml(k)}</td>
                <td class="py-1 pr-3 text-white">${w}</td>
                <td class="py-1 text-slate-400">${escapeHtml(meta.component_descriptions[k] || "")}</td>
              </tr>`).join("")}</tbody>
          </table>
          <p class="mt-2 text-amber-300/80">${escapeHtml(meta.note)}</p>`;
    } catch (e) {
        box.textContent = `Failed to load methodology: ${e.message}`;
    }
}

async function init() {
    try {
        const drugs = await apiGet("/drugs");
        const withCandidates = await apiGet("/interactions/candidates");
        const idsWithCandidates = new Set(withCandidates.map((it) => it.drug_id));
        const drugSelect = document.getElementById("flt-drug");
        drugSelect.innerHTML = '<option value="">All drugs</option>' +
            drugs.filter((d) => idsWithCandidates.has(d.drug_id))
                .map((d) => `<option value="${escapeHtml(d.drug_id)}">${escapeHtml(d.name)}</option>`).join("");
        await Promise.all([render(), loadMethodology()]);
    } catch (e) {
        document.getElementById("cand-status").textContent = `Failed to load: ${e.message}`;
    }
}

document.getElementById("flt-drug").addEventListener("change", () => {
    state.drug_id = document.getElementById("flt-drug").value;
    render();
});
document.getElementById("flt-sort").addEventListener("change", () => {
    state.sort = document.getElementById("flt-sort").value;
    render();
});
document.getElementById("flt-low-resource").addEventListener("change", () => {
    state.low_resource_relevance = document.getElementById("flt-low-resource").checked;
    render();
});
let debounce;
document.getElementById("flt-tissue").addEventListener("input", () => {
    clearTimeout(debounce);
    debounce = setTimeout(() => {
        state.tissue_origin = document.getElementById("flt-tissue").value.trim();
        render();
    }, 250);
});

init();
