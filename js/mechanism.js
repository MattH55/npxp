/* Mechanism-Based Hypotheses page. Real, cited lookups (SynLethDB +
   DGIdb), not a learned score -- see mechanism-hypotheses.html's own
   banner and pipeline/scripts/build_mechanism_hypotheses.py. */
renderNav("Mechanism Hypotheses");

const TIER_LABEL = { 3: "CRISPR/CRISPRi (direct screen)", 2: "Curated / high-throughput", 1: "Computational prediction", 0: "Text mining" };
const TIER_CLASS = { 3: "bg-emerald-600", 2: "bg-sky-600", 1: "bg-amber-600", 0: "bg-slate-600" };
const CATEGORY_LABEL = {
    induced_hr_deficiency: "Induced HR deficiency → synthetic lethality",
    pathway_suppression: "Pathway suppression → same-axis drug class",
    immune_mobilization: "Immune mobilisation → checkpoint inhibition",
};

let ALL_ROWS = [];
let MODIFIER_ORDER = [];
let ACTIVE_MODIFIER = null;

function citationString(c) {
    // {pmid, first_author, journal, year} -> "pmid:XXXXX" for citationChips,
    // with the readable label built separately since citationChips() only
    // shows "PMID XXXXX" by default.
    return `pmid:${c.pmid}`;
}

function citationLine(c, note) {
    const [href] = citationHref(citationString(c));
    const label = `${escapeHtml(c.first_author)} et al., <i>${escapeHtml(c.journal)}</i> ${c.year}`;
    const link = href ? `<a class="underline hover:text-white" target="_blank" rel="noopener" href="${href}">${label} (PMID ${c.pmid})</a>` : label;
    return `<div class="text-xs text-slate-400">${link}${note ? ` — ${escapeHtml(note)}` : ""}</div>`;
}

function groupByDrug(rows) {
    const byDrug = new Map();
    for (const r of rows) {
        if (!byDrug.has(r.drug_name)) {
            byDrug.set(r.drug_name, {
                drug_name: r.drug_name, targets: new Set(), best_tier: -1,
                n_sources: 0, sl_citations: [], has_direct: false, direct_citations: [],
            });
        }
        const g = byDrug.get(r.drug_name);
        g.targets.add(r.drug_target_gene);
        g.n_sources = Math.max(g.n_sources, r.drug_n_sources);
        const tier = r.sl_evidence_tier ?? 2;
        if (tier > g.best_tier) g.best_tier = tier;
        if (r.sl_pubmed_id) g.sl_citations.push({ gene: r.drug_target_gene, pmid: r.sl_pubmed_id, evidence: r.sl_evidence_type });
        if (r.has_direct_experimental_support) {
            g.has_direct = true;
            g.direct_citations.push(...r.direct_experimental_support);
        }
    }
    return [...byDrug.values()].sort((a, b) => (b.has_direct - a.has_direct) || (b.best_tier - a.best_tier) || (b.n_sources - a.n_sources));
}

function renderTabs() {
    const el = document.getElementById("tabs");
    el.innerHTML = MODIFIER_ORDER.map((m) => {
        const active = m.modifier_id === ACTIVE_MODIFIER;
        const cls = active ? "bg-indigo-600 text-white" : "bg-slate-800 text-slate-300 hover:text-white";
        return `<button data-modifier="${m.modifier_id}" class="px-3 py-1.5 text-sm rounded-lg font-semibold ${cls}">${escapeHtml(m.modifier_label)}</button>`;
    }).join("");
    el.querySelectorAll("button").forEach((b) => b.addEventListener("click", () => {
        ACTIVE_MODIFIER = b.dataset.modifier;
        renderTabs();
        renderContent();
    }));
}

function renderContent() {
    const rows = ALL_ROWS.filter((r) => r.modifier_id === ACTIVE_MODIFIER);
    if (!rows.length) {
        document.getElementById("content").innerHTML = `<div class="text-slate-500 text-sm">No real hypotheses resolved for this modifier.</div>`;
        return;
    }
    const first = rows[0];
    const grouped = groupByDrug(rows);

    const header = `
      <div class="bg-slate-800/40 border border-slate-700 rounded-lg p-4 mb-4">
        <div class="text-xs uppercase tracking-wide text-indigo-400 font-semibold mb-1">${escapeHtml(CATEGORY_LABEL[first.category] || first.category)}</div>
        <p class="text-sm text-slate-300 mb-2">${escapeHtml(first.mechanism)}</p>
        ${citationLine(first.modifier_citation)}
        ${(first.supporting_citations || []).map((c) => citationLine(c, c.note)).join("")}
      </div>`;

    const rowsHtml = grouped.map((g) => `
      <tr class="bg-slate-800/40">
        <td class="py-2 pl-2 rounded-l-lg">
          <div class="text-white font-semibold">${escapeHtml(g.drug_name)}</div>
          ${g.has_direct ? `<span class="bg-amber-500 text-slate-900 text-xs font-bold px-2 py-0.5 rounded-full"><i class="fa-solid fa-star mr-1"></i>Direct experimental support</span>` : ""}
        </td>
        <td class="text-slate-300 text-sm">${[...g.targets].map(escapeHtml).join(", ")}</td>
        <td><span class="${TIER_CLASS[g.best_tier] ?? "bg-slate-600"} text-white text-xs font-semibold px-2 py-0.5 rounded-full">${escapeHtml(TIER_LABEL[g.best_tier] ?? "Literature-direct")}</span></td>
        <td class="text-slate-300 text-sm">${g.n_sources} source${g.n_sources === 1 ? "" : "s"}</td>
        <td class="text-xs text-slate-400">
          ${g.direct_citations.map((c) => citationLine(c, c.note)).join("")}
          ${g.sl_citations.slice(0, 3).map((c) => c.pmid ? `<div>SL: ${escapeHtml(c.gene)} (PMID <a class="underline hover:text-white" target="_blank" href="https://pubmed.ncbi.nlm.nih.gov/${c.pmid}/">${c.pmid}</a>, ${escapeHtml(c.evidence)})</div>` : "").join("")}
        </td>
      </tr>`).join("");

    document.getElementById("content").innerHTML = header + `
      <table class="w-full text-sm border-separate" style="border-spacing:0 4px">
        <thead><tr class="text-left text-xs text-slate-400">
          <th class="py-1 pl-2">Drug</th>
          <th class="py-1">Target gene(s)</th>
          <th class="py-1">Best evidence</th>
          <th class="py-1">DGIdb corroboration</th>
          <th class="py-1">Citations</th>
        </tr></thead>
        <tbody>${rowsHtml}</tbody>
      </table>`;
}

async function main() {
    const statusEl = document.getElementById("status");
    statusEl.textContent = "Loading...";
    try {
        const data = await apiGet("/mechanism-hypotheses");
        ALL_ROWS = data.rows;
        const seen = new Set();
        MODIFIER_ORDER = [];
        for (const r of ALL_ROWS) {
            if (seen.has(r.modifier_id)) continue;
            seen.add(r.modifier_id);
            MODIFIER_ORDER.push({ modifier_id: r.modifier_id, modifier_label: r.modifier_label });
        }
        ACTIVE_MODIFIER = MODIFIER_ORDER[0]?.modifier_id;
        statusEl.textContent = `${data.n_rows} real, cited hypothesis rows across ${MODIFIER_ORDER.length} modalities. Built ${new Date(data.built).toLocaleString()}.`;
        renderTabs();
        renderContent();
    } catch (e) {
        statusEl.textContent = `Failed to load: ${e.message}`;
    }
}

main();
