/* Explorer: drug selector + cell-line × modifier interaction heatmap. */
renderNav("Explorer");

const TYPE_CLASS = { synergistic: "syn", additive: "add", antagonistic: "ant", unknown: "unk" };
const TIER_CLASS = {
    tier_1_direct: "t1",
    tier_2_inferred: "t2",
    tier_2b_model_predicted: "t2b",
    tier_3_mechanism_only: "t3",
};
const MOD_TYPE_LABEL = {
    dietary_metabolic: "Dietary/metabolic",
    thermal: "Thermal",
    hypoxic: "Hypoxic",
    mechanical_radiative: "Mechanical/radiative",
    other: "Other",
};

function cellLabel(cell) {
    if (cell.combined_effect_metric !== null && cell.combined_effect_metric !== undefined) {
        return Number(cell.combined_effect_metric).toFixed(2);
    }
    return cell.interaction_type === "unknown" ? "?" : cell.interaction_type.slice(0, 3).toUpperCase();
}

function modifierColHeader(m) {
    const params = m.protocol_parameters || {};
    const bits = [];
    if (params.temperature_C) bits.push(`${params.temperature_C} °C`);
    if (params.concentration_mM) bits.push(`${params.concentration_mM} mM`);
    if (params.o2_percent) bits.push(`${params.o2_percent}% O₂`);
    if (params.duration_hr) bits.push(`${params.duration_hr} h`);
    if (params.duration_min) bits.push(`${params.duration_min} min`);
    return `<div class="text-white text-xs font-semibold leading-tight">${escapeHtml(m.agent)}</div>
            <div class="text-slate-500 text-[0.65rem]">${escapeHtml(MOD_TYPE_LABEL[m.modifier_type] || m.modifier_type)}</div>
            <div class="text-slate-400 text-[0.65rem]">${escapeHtml(bits.join(" · "))}</div>`;
}

function renderMatrix(data) {
    const status = document.getElementById("explorer-status");
    if (!data.cells.length) {
        status.textContent = "No curated interactions for this selection yet.";
        document.getElementById("heatmap").innerHTML = "";
        return;
    }
    const drugName = data.drug ? data.drug.name : (data.drug_class || "");
    status.textContent = `${data.cells.length} interaction claim(s) · ${data.cell_lines.length} cell line(s) · ${data.modifiers.length} modifier(s)`;

    const cellMap = new Map(data.cells.map((c) => [`${c.cell_line_id}|${c.modifier_id}`, c]));
    let html = `<table class="border-separate" style="border-spacing:6px"><thead><tr>
        <th class="text-left align-bottom pb-2 pr-2 text-slate-400 text-xs font-medium">
          ${escapeHtml(drugName)}<div class="text-slate-600">cell line ↓ · modifier →</div></th>`;
    for (const m of data.modifiers) {
        html += `<th class="align-bottom pb-2" style="min-width:110px">${modifierColHeader(m)}</th>`;
    }
    html += `</tr></thead><tbody>`;

    for (const cl of data.cell_lines) {
        html += `<tr>
          <td class="pr-2">
            <a href="cell-lines.html#${encodeURIComponent(cl.cell_line_id)}" class="text-white text-sm font-semibold hover:text-emerald-400">${escapeHtml(cl.name)}</a>
            <div class="text-slate-500 text-[0.65rem]">${escapeHtml(cl.tissue_origin)}</div>
          </td>`;
        for (const m of data.modifiers) {
            const cell = cellMap.get(`${cl.cell_line_id}|${m.modifier_id}`);
            if (!cell) {
                html += `<td><div class="hm-cell none">—</div></td>`;
                continue;
            }
            const tc = TYPE_CLASS[cell.interaction_type] || "unk";
            const tierCls = TIER_CLASS[cell.evidence_tier] || "t3";
            const tier = TIER_SHORT[cell.evidence_tier] || "?";
            const title = `${cl.name} × ${m.agent} + ${drugName} — ${cell.interaction_type} (${tier})`;
            html += `<td>
              <div class="hm-cell ${tc} ${tierCls}" title="${escapeHtml(title)}"
                   onclick="location.href='interaction.html?id=${cell.interaction_id}'">
                <span class="font-bold">${escapeHtml(cellLabel(cell))}</span>
                <span class="tier">${tier}</span>
              </div></td>`;
        }
        html += `</tr>`;
    }
    html += `</tbody></table>`;
    document.getElementById("heatmap").innerHTML = html;
}

async function loadMatrix(drugId) {
    const status = document.getElementById("explorer-status");
    status.textContent = "Loading…";
    try {
        const data = await apiGet("/explorer/matrix", { drug_id: drugId });
        renderMatrix(data);
    } catch (e) {
        status.textContent = `Failed to load matrix: ${e.message}`;
        document.getElementById("heatmap").innerHTML = "";
    }
}

(async () => {
    const select = document.getElementById("drug-select");
    try {
        const drugs = await apiGet("/drugs");
        const withData = drugs.filter((d) => d.interaction_count > 0);
        select.innerHTML =
            withData.map((d) =>
                `<option value="${escapeHtml(d.drug_id)}">${escapeHtml(d.name)}${d.drug_class ? ` (${escapeHtml(d.drug_class)})` : ""} — ${d.interaction_count}</option>`
            ).join("");
        const params = new URLSearchParams(location.search);
        const requested = params.get("drug_id");
        const initial = withData.some((d) => d.drug_id === requested) ? requested
            : (withData[0] ? withData[0].drug_id : null);
        if (!initial) {
            document.getElementById("explorer-status").textContent = "No interactions curated yet.";
            return;
        }
        select.value = initial;
        select.addEventListener("change", () => loadMatrix(select.value));
        await loadMatrix(initial);
    } catch (e) {
        document.getElementById("explorer-status").textContent = `Failed to load drugs: ${e.message}`;
    }
})();
