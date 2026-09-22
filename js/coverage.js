/* Coverage Gap Matrix (build spec v6): drugs x modifier_type, shaded by
   gap_score. Rendered as a real table, not a canvas heatmap, so cell values
   stay screen-reader/copy-paste accessible. */
renderNav("Coverage Gaps");

const MOD_TYPE_LABEL = {
    dietary_metabolic: "Dietary/metabolic",
    thermal: "Thermal",
    hypoxic: "Hypoxic",
    mechanical_radiative: "Mechanical/radiative",
    other: "Other",
};

function shade(gapScore, maxScore) {
    if (maxScore <= 0) return "rgba(100,116,139,0.15)";
    const t = Math.min(gapScore / maxScore, 1);
    // amber ramp: low gap -> subtle, high gap -> saturated.
    const alpha = 0.12 + 0.68 * t;
    return `rgba(245,158,11,${alpha.toFixed(2)})`;
}

async function render(tissueOrigin) {
    const status = document.getElementById("cov-status");
    status.textContent = "Loading…";
    try {
        const rows = await apiGet("/coverage-gaps", tissueOrigin ? { tumor_type: tissueOrigin } : {});
        if (!rows.length) {
            status.textContent = "No drugs with resistance-mechanism annotation yet.";
            document.getElementById("cov-matrix").innerHTML = "";
            return;
        }
        const maxScore = Math.max(...rows.map((r) => r.gap_score));
        const drugIds = [...new Set(rows.map((r) => r.drug_id))];
        const drugNames = Object.fromEntries(rows.map((r) => [r.drug_id, r.drug_name]));
        const modTypes = [...new Set(rows.map((r) => r.modifier_type))];
        const byKey = new Map(rows.map((r) => [`${r.drug_id}|${r.modifier_type}`, r]));

        status.textContent = `${drugIds.length} drug(s) x ${modTypes.length} modifier type(s) — max gap_score ${maxScore.toFixed(2)}`;

        let html = `<table class="border-separate" style="border-spacing:4px"><thead><tr>
          <th class="text-left align-bottom pb-2 pr-2 text-slate-400 text-xs font-medium">drug ↓ · modifier type →</th>`;
        for (const mt of modTypes) {
            html += `<th class="align-bottom pb-2 text-slate-300 text-xs font-semibold" style="min-width:120px">${escapeHtml(MOD_TYPE_LABEL[mt] || mt)}</th>`;
        }
        html += `</tr></thead><tbody>`;
        for (const drugId of drugIds) {
            html += `<tr><td class="pr-2"><a href="drugs.html#${encodeURIComponent(drugId)}" class="text-white text-sm font-semibold hover:text-emerald-400">${escapeHtml(drugNames[drugId])}</a></td>`;
            for (const mt of modTypes) {
                const cell = byKey.get(`${drugId}|${mt}`);
                if (!cell) { html += `<td></td>`; continue; }
                const bg = shade(cell.gap_score, maxScore);
                const title = `${drugNames[drugId]} x ${MOD_TYPE_LABEL[mt] || mt}: gap_score ${cell.gap_score} (${cell.mechanism_count} mechanism(s), ${cell.interaction_count} interaction(s))`;
                html += `<td>
                  <div class="rounded-lg p-2 text-center text-xs" style="background:${bg}" title="${escapeHtml(title)}">
                    <div class="text-white font-bold">${cell.gap_score}</div>
                    <div class="text-slate-300">${cell.interaction_count} tested</div>
                  </div>
                </td>`;
            }
            html += `</tr>`;
        }
        html += `</tbody></table>`;
        document.getElementById("cov-matrix").innerHTML = html;
    } catch (e) {
        status.textContent = `Failed to load coverage gaps: ${e.message}`;
    }
}

let debounce;
document.getElementById("flt-tissue").addEventListener("input", () => {
    clearTimeout(debounce);
    debounce = setTimeout(() => render(document.getElementById("flt-tissue").value.trim()), 250);
});

render();
