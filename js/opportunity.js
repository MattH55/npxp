/* Opportunity Ranking (reconciling the "NPxP Interaction Predictor" build
   spec's section 6): cell-line and cancer-type views over
   /api/opportunity/cell-lines and /api/opportunity/cancer-types. A
   genuinely different question/formula from candidates.js's row-level
   priority_score worklist -- see opportunity.html's own explanation and
   synlethality/opportunity_scoring.py. */
renderNav("Opportunity");

let cellLineRows = [];
let cancerTypeRows = [];
let cellLineSort = { key: "opportunity_score", dir: -1 };
let cancerTypeSort = { key: "mean_opportunity_score", dir: -1 };

function sortRows(rows, sort) {
    const copy = [...rows];
    copy.sort((a, b) => {
        const av = a[sort.key], bv = b[sort.key];
        if (typeof av === "string") return sort.dir * av.localeCompare(bv);
        return sort.dir * ((av ?? -Infinity) - (bv ?? -Infinity));
    });
    return copy;
}

function sortHeader(label, key, currentSort, onClick) {
    const active = currentSort.key === key;
    const arrow = active ? (currentSort.dir === -1 ? "▼" : "▲") : "";
    return `<th class="py-1 pr-3 cursor-pointer select-none text-slate-400 hover:text-white ${active ? "text-white" : ""}"
              data-sort-key="${key}">${escapeHtml(label)} ${arrow}</th>`;
}

function renderCellLines() {
    const rows = sortRows(cellLineRows, cellLineSort);
    const html = `
      <table class="w-full text-sm border-separate" style="border-spacing:0 4px">
        <thead><tr class="text-left text-xs">
          ${sortHeader("Cell line", "cell_line_name", cellLineSort)}
          ${sortHeader("Tissue", "tissue_origin", cellLineSort)}
          ${sortHeader("Cancer type", "oncotree_primary_disease", cellLineSort)}
          ${sortHeader("Opportunity", "opportunity_score", cellLineSort)}
          ${sortHeader("Predicted uplift", "predicted_uplift", cellLineSort)}
          ${sortHeader("Gap weight", "gap_weight", cellLineSort)}
          ${sortHeader("Tier 1 rows", "tier_1_direct_count", cellLineSort)}
          ${sortHeader("All real rows", "n_real_interactions", cellLineSort)}
        </tr></thead>
        <tbody>
          ${rows.map((r) => `
            <tr class="bg-slate-800/40">
              <td class="py-2 pl-2 rounded-l-lg"><a class="text-white font-semibold hover:text-emerald-400" href="cell-lines.html#${encodeURIComponent(r.cell_line_id)}">${escapeHtml(r.cell_line_name)}</a></td>
              <td class="text-slate-300">${escapeHtml(r.tissue_origin)}</td>
              <td class="text-slate-400">${escapeHtml(r.oncotree_primary_disease || "—")}</td>
              <td class="text-indigo-300 font-bold">${r.opportunity_score.toFixed(3)}</td>
              <td class="text-slate-300">${r.predicted_uplift.toFixed(3)}</td>
              <td class="text-slate-300">${r.gap_weight.toFixed(3)}</td>
              <td class="text-slate-300">${r.tier_1_direct_count}</td>
              <td class="text-slate-300 rounded-r-lg pr-2">${r.n_real_interactions} / ${r.n_total_interactions} total</td>
            </tr>`).join("")}
        </tbody>
      </table>`;
    document.getElementById("opp-cell-lines").innerHTML = html;
    document.querySelectorAll("#opp-cell-lines [data-sort-key]").forEach((th) => {
        th.addEventListener("click", () => {
            const key = th.dataset.sortKey;
            cellLineSort.dir = cellLineSort.key === key ? -cellLineSort.dir : -1;
            cellLineSort.key = key;
            renderCellLines();
        });
    });
}

function renderCancerTypes() {
    const rows = sortRows(cancerTypeRows, cancerTypeSort);
    const html = `
      <table class="w-full text-sm border-separate" style="border-spacing:0 4px">
        <thead><tr class="text-left text-xs">
          ${sortHeader("Cancer type", "oncotree_primary_disease", cancerTypeSort)}
          ${sortHeader("Mean opportunity", "mean_opportunity_score", cancerTypeSort)}
          ${sortHeader("Max opportunity", "max_opportunity_score", cancerTypeSort)}
          ${sortHeader("Cell lines", "n_cell_lines", cancerTypeSort)}
          <th class="py-1 pr-3 text-slate-400">Confidence</th>
        </tr></thead>
        <tbody>
          ${rows.map((r) => `
            <tr class="bg-slate-800/40">
              <td class="py-2 pl-2 rounded-l-lg text-white font-semibold">${escapeHtml(r.oncotree_primary_disease)}</td>
              <td class="text-indigo-300 font-bold">${r.mean_opportunity_score.toFixed(3)}</td>
              <td class="text-slate-300">${r.max_opportunity_score.toFixed(3)}</td>
              <td class="text-slate-300" title="${r.cell_line_ids.map(escapeHtml).join(', ')}">${r.n_cell_lines}</td>
              <td class="rounded-r-lg pr-2">${r.low_confidence
                ? '<span class="text-amber-400"><i class="fa-solid fa-triangle-exclamation mr-1"></i>Low (1 cell line)</span>'
                : '<span class="text-slate-400">OK</span>'}</td>
            </tr>`).join("")}
        </tbody>
      </table>`;
    document.getElementById("opp-cancer-types").innerHTML = html;
    document.querySelectorAll("#opp-cancer-types [data-sort-key]").forEach((th) => {
        th.addEventListener("click", () => {
            const key = th.dataset.sortKey;
            cancerTypeSort.dir = cancerTypeSort.key === key ? -cancerTypeSort.dir : -1;
            cancerTypeSort.key = key;
            renderCancerTypes();
        });
    });
}

function showTab(which) {
    const cellLinesTab = document.getElementById("tab-cell-lines");
    const cancerTypesTab = document.getElementById("tab-cancer-types");
    const cellLinesPanel = document.getElementById("opp-cell-lines");
    const cancerTypesPanel = document.getElementById("opp-cancer-types");
    const active = "bg-indigo-600 text-white font-semibold";
    const inactive = "bg-slate-800 text-slate-300";
    if (which === "cell-lines") {
        cellLinesTab.className = `px-3 py-1.5 text-sm rounded-lg ${active}`;
        cancerTypesTab.className = `px-3 py-1.5 text-sm rounded-lg ${inactive}`;
        cellLinesPanel.classList.remove("hidden");
        cancerTypesPanel.classList.add("hidden");
    } else {
        cancerTypesTab.className = `px-3 py-1.5 text-sm rounded-lg ${active}`;
        cellLinesTab.className = `px-3 py-1.5 text-sm rounded-lg ${inactive}`;
        cancerTypesPanel.classList.remove("hidden");
        cellLinesPanel.classList.add("hidden");
    }
}

document.getElementById("tab-cell-lines").addEventListener("click", () => showTab("cell-lines"));
document.getElementById("tab-cancer-types").addEventListener("click", () => showTab("cancer-types"));

(async () => {
    const status = document.getElementById("opp-status");
    status.textContent = "Loading…";
    try {
        [cellLineRows, cancerTypeRows] = await Promise.all([
            apiGet("/opportunity/cell-lines"),
            apiGet("/opportunity/cancer-types"),
        ]);
        status.textContent = `${cellLineRows.length} cell line(s) across ${cancerTypeRows.length} cancer type grouping(s).`;
        renderCellLines();
        renderCancerTypes();
    } catch (e) {
        status.textContent = `Failed to load opportunity scores: ${e.message}`;
    }
})();
