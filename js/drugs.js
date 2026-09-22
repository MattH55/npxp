/* Drug library: filterable list + detail panel with cross-refs and known
   interactions/responses. Client-side filtering here matches the pattern in
   modifiers.js; the live backend does its own filtering server-side via
   /api/drugs (list_drugs has no query params today, so filters are applied
   client-side here in both live and static modes for consistency). */
renderNav("Drugs");

const STATUS_LABEL = {
    approved: "Approved",
    investigational: "Investigational",
    preclinical: "Preclinical",
    withdrawn: "Withdrawn",
};
const STATUS_STYLE = {
    approved: "bg-emerald-600",
    investigational: "bg-amber-600",
    preclinical: "bg-slate-600",
    withdrawn: "bg-rose-700",
};

function statusBadge(status) {
    if (!status) return '<span class="text-slate-500 text-xs">unspecified</span>';
    const cls = STATUS_STYLE[status] || "bg-slate-600";
    return `<span class="${cls} text-white text-xs font-semibold px-2 py-0.5 rounded-full">${escapeHtml(STATUS_LABEL[status] || status)}</span>`;
}

function crossRefChips(d) {
    const chips = [];
    if (d.pubchem_cid) chips.push(`<a class="cite-chip" target="_blank" rel="noopener" href="https://pubchem.ncbi.nlm.nih.gov/compound/${encodeURIComponent(d.pubchem_cid)}">PubChem ${escapeHtml(d.pubchem_cid)} <i class="fa-solid fa-arrow-up-right-from-square text-[0.6rem]"></i></a>`);
    if (d.chembl_id) chips.push(`<a class="cite-chip" target="_blank" rel="noopener" href="https://www.ebi.ac.uk/chembl/compound_report_card/${encodeURIComponent(d.chembl_id)}/">ChEMBL ${escapeHtml(d.chembl_id)} <i class="fa-solid fa-arrow-up-right-from-square text-[0.6rem]"></i></a>`);
    if (d.drugbank_id) chips.push(`<a class="cite-chip" target="_blank" rel="noopener" href="https://go.drugbank.com/drugs/${encodeURIComponent(d.drugbank_id)}">DrugBank ${escapeHtml(d.drugbank_id)} <i class="fa-solid fa-arrow-up-right-from-square text-[0.6rem]"></i></a>`);
    return chips.join("") || '<span class="text-slate-500 text-xs">No verified cross-references yet</span>';
}

let allDrugs = [];
const state = { clinical_status: "", q: "" };

function applyFilters(rows) {
    const q = state.q.toLowerCase();
    return rows.filter((d) =>
        (!state.clinical_status || d.clinical_status === state.clinical_status) &&
        (!q || [d.name, d.drug_class, d.target, d.mechanism_of_action].some(
            (f) => (f || "").toLowerCase().includes(q))));
}

async function loadList() {
    const status = document.getElementById("drug-status");
    status.textContent = "Loading…";
    try {
        if (!allDrugs.length) allDrugs = await apiGet("/drugs");
        const rows = applyFilters(allDrugs);
        status.textContent = `${rows.length} of ${allDrugs.length} drug(s)`;
        document.getElementById("drug-list").innerHTML = rows.map((d) => `
          <div class="row-link bg-slate-800/60 border border-slate-700 rounded-xl p-4" data-id="${escapeHtml(d.drug_id)}">
            <div class="flex items-start justify-between gap-2">
              <div class="text-white font-semibold">${escapeHtml(d.name)}</div>
              ${statusBadge(d.clinical_status)}
            </div>
            <div class="text-slate-400 text-sm mt-1">${escapeHtml(d.drug_class || "—")}</div>
            <div class="text-slate-500 text-xs mt-1">${escapeHtml(d.target || "")}</div>
            <div class="flex items-center justify-between mt-2 text-xs text-slate-500">
              <span>${d.interaction_count} interaction(s) · ${d.response_count} response(s)</span>
              ${d.induced_expression_signature_ref ? '<span class="cite-chip">LINCS-backed</span>' : ""}
            </div>
          </div>`).join("") || '<div class="text-slate-500 text-sm">No drugs match.</div>';
        document.querySelectorAll("#drug-list .row-link").forEach((el) =>
            el.addEventListener("click", () => {
                location.hash = el.dataset.id;
                loadDetail(el.dataset.id);
            }));
    } catch (e) {
        status.textContent = `Failed to load drugs: ${e.message}`;
    }
}

async function loadDetail(drugId) {
    const box = document.getElementById("drug-detail");
    box.innerHTML = `<div class="text-slate-400 text-sm">Loading ${escapeHtml(drugId)}…</div>`;
    try {
        const d = await apiGet(`/drugs/${encodeURIComponent(drugId)}`);
        const intRows = (d.interactions || []).map((it) => `
          <tr class="row-link border-b border-slate-700/60" onclick="location.href='interaction.html?id=${it.id}'">
            <td class="py-2 pr-3 text-white">${escapeHtml(it.modifier_id)}${it.is_bridging_candidate ? ' <i class="fa-solid fa-flask-vial text-[0.6rem] text-violet-400" title="Mechanistic Bridging candidate"></i>' : ""}</td>
            <td class="py-2 pr-3 text-slate-300">${escapeHtml(it.cell_line_name)}</td>
            <td class="py-2 pr-3">${typeBadge(it.interaction_type)}</td>
            <td class="py-2 pr-3">${tierBadge(it.evidence_tier)}</td>
            <td class="py-2">${citationChips(it.source_study)}</td>
          </tr>`).join("");
        const respRows = (d.drug_responses || []).map((r) => `
          <tr class="border-b border-slate-700/60">
            <td class="py-2 pr-3 text-slate-300">${escapeHtml(r.cell_line_id)}</td>
            <td class="py-2 pr-3 text-slate-300">${escapeHtml(r.metric_type)}</td>
            <td class="py-2 pr-3 text-slate-300">${r.viability_metric}</td>
            <td class="py-2 text-slate-400">${escapeHtml(r.source)}</td>
          </tr>`).join("");
        box.innerHTML = `
          <div class="bg-slate-800/60 border border-slate-700 rounded-xl p-6">
            <div class="flex items-start justify-between gap-2">
              <h2 class="text-xl font-bold text-white mb-1">${escapeHtml(d.name)}
                <span class="text-slate-500 text-sm font-normal ml-2">${escapeHtml(d.drug_id)}</span></h2>
              ${statusBadge(d.clinical_status)}
            </div>
            <div class="text-slate-400 text-sm mb-1">${escapeHtml(d.drug_class || "—")}</div>
            <div class="text-slate-400 text-sm mb-3">${escapeHtml(d.mechanism_of_action || d.target || "")}</div>
            <div class="flex flex-wrap gap-2 mb-2">${crossRefChips(d)}</div>
            <div class="text-xs text-slate-400 mb-4">Cost/accessibility tier: <span class="cite-chip">${escapeHtml(d.cost_accessibility_tier || "unassessed")}</span></div>
            ${(d.synonyms || []).length ? `<div class="text-slate-500 text-xs mb-4">Also known as: ${d.synonyms.map(escapeHtml).join(", ")}</div>` : ""}
            ${d.induced_expression_signature_ref
                ? `<div class="text-xs mb-4"><span class="cite-chip">LINCS signature: ${escapeHtml(d.induced_expression_signature_ref)}</span></div>`
                : `<div class="text-slate-500 text-xs mb-4">No LINCS L1000 signature yet — Step 2 prediction cannot use this drug's induced-expression profile until the ingestion pipeline is enabled.</div>`}
            ${(d.resistance_mechanisms || []).length ? `
              <h3 class="text-white font-semibold mt-2"><i class="fa-solid fa-shield-halved mr-2 text-emerald-400"></i>Resistance mechanisms</h3>
              <div class="space-y-2 mt-2 mb-4">${d.resistance_mechanisms.map((rm) => `
                <div class="bg-slate-900/60 border border-slate-700 rounded-lg p-3 text-sm">
                  <div class="text-white font-semibold">${escapeHtml(rm.pathway_or_gene)}</div>
                  <div class="text-slate-300 mt-1">${escapeHtml(rm.mechanism_description)}</div>
                  <div class="mt-2">${citationChips([rm.source])}</div>
                </div>`).join("")}</div>` : ""}
            <h3 class="text-white font-semibold mt-2"><i class="fa-solid fa-circle-nodes mr-2 text-emerald-400"></i>Modifier interactions</h3>
            ${intRows ? `<table class="w-full text-sm mt-2">
              <thead><tr class="text-left text-slate-400 border-b border-slate-700">
                <th class="py-2 pr-3">Modifier</th><th class="py-2 pr-3">Cell line</th>
                <th class="py-2 pr-3">Type</th><th class="py-2 pr-3">Tier</th><th class="py-2">Citations</th>
              </tr></thead><tbody>${intRows}</tbody></table>` : '<div class="text-slate-500 text-sm">None recorded.</div>'}
            <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-vial mr-2 text-emerald-400"></i>Monotherapy responses (DepMap/PRISM, GDSC, CTRP)</h3>
            ${respRows ? `<table class="w-full text-sm mt-2">
              <thead><tr class="text-left text-slate-400 border-b border-slate-700">
                <th class="py-2 pr-3">Cell line</th><th class="py-2 pr-3">Metric</th>
                <th class="py-2 pr-3">Value</th><th class="py-2">Source</th>
              </tr></thead><tbody>${respRows}</tbody></table>` : '<div class="text-slate-500 text-sm">None ingested yet (ingestion pipelines are structured stubs pending source data review).</div>'}
          </div>`;
        box.scrollIntoView({ behavior: "smooth", block: "start" });
    } catch (e) {
        box.innerHTML = `<div class="text-rose-400 text-sm">${escapeHtml(e.message)}</div>`;
    }
}

let debounce;
document.getElementById("flt-status").addEventListener("change", () => {
    state.clinical_status = document.getElementById("flt-status").value;
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
