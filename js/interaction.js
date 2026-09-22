/* Interaction detail page: drug + modifier + cell line, mechanism, citations, tier. */
renderNav("Explorer");

function mechanismBlock(mech) {
    if (!mech) {
        return `<div class="text-slate-500 text-sm">
          No stress-signature mechanism is linked to this claim
          (Tier 1 direct evidence does not require a mechanistic anchor).</div>`;
    }
    return `<table class="w-full text-sm mt-2">
      <tbody>
        <tr class="border-b border-slate-700/60">
          <td class="py-2 pr-3 kv-key w-40">Panel</td><td class="py-2 text-slate-200">${escapeHtml(mech.signature_panel)}</td></tr>
        <tr class="border-b border-slate-700/60">
          <td class="py-2 pr-3 kv-key">Directionality</td><td class="py-2 text-slate-200">${escapeHtml(mech.directionality)}</td></tr>
        <tr class="border-b border-slate-700/60">
          <td class="py-2 pr-3 kv-key">Score</td><td class="py-2 text-slate-200">${mech.score === null ? '<span class="text-slate-500">n/a (see notes)</span>' : escapeHtml(mech.score)}</td></tr>
        <tr class="border-b border-slate-700/60">
          <td class="py-2 pr-3 kv-key">Evidence ref</td><td class="py-2 text-slate-200">${escapeHtml(mech.raw_deg_evidence_ref)}</td></tr>
        <tr><td class="py-2 pr-3 kv-key">Notes</td><td class="py-2 text-slate-300">${escapeHtml(mech.notes)}</td></tr>
      </tbody></table>`;
}

function relatedBlock(related) {
    if (!related.length) return "";
    return `<h3 class="text-white font-semibold mt-8 mb-2">
        <i class="fa-solid fa-circle-nodes mr-2 text-emerald-400"></i>Other interactions in this cell line</h3>
      <table class="w-full text-sm">
        <thead><tr class="text-left text-slate-400 border-b border-slate-700">
          <th class="py-2 pr-3">Drug</th><th class="py-2 pr-3">Modifier</th>
          <th class="py-2 pr-3">Type</th><th class="py-2 pr-3">Tier</th>
        </tr></thead>
        <tbody>${related.map((r) => `
          <tr class="row-link border-b border-slate-700/60" onclick="location.href='interaction.html?id=${r.id}'">
            <td class="py-2 pr-3 text-white">${escapeHtml(r.drug_name)}</td>
            <td class="py-2 pr-3 text-slate-300">${escapeHtml(r.modifier_agent)}</td>
            <td class="py-2 pr-3">${typeBadge(r.interaction_type)}</td>
            <td class="py-2 pr-3">${tierBadge(r.evidence_tier)}</td>
          </tr>`).join("")}
        </tbody></table>`;
}

async function load() {
    const box = document.getElementById("interaction-detail");
    const id = new URLSearchParams(location.search).get("id");
    if (!id) {
        box.innerHTML = '<div class="text-rose-400 text-sm">Missing ?id= parameter.</div>';
        return;
    }
    box.innerHTML = '<div class="text-slate-400 text-sm">Loading…</div>';
    try {
        const it = await apiGet(`/interactions/${encodeURIComponent(id)}`);
        const tierDef = it.evidence_tier_definition || {};
        box.innerHTML = `
          <div class="bg-slate-800/60 border border-slate-700 rounded-xl p-6">
            <div class="flex flex-wrap items-center gap-2 mb-3">
              ${typeBadge(it.interaction_type)} ${tierBadge(it.evidence_tier)}
              <span class="text-slate-500 text-xs ml-1">${escapeHtml(tierDef.label || "")}</span>
            </div>
            <h1 class="text-2xl font-extrabold text-white leading-snug">
              ${escapeHtml(it.drug_name)}
              <span class="text-slate-500 font-normal">+</span> ${escapeHtml(it.modifier_agent)}
              <span class="text-slate-500 font-normal">in</span> ${escapeHtml(it.cell_line_name)}
            </h1>
            <div class="text-slate-400 text-sm mt-1">${escapeHtml(it.tissue_origin || "")} ·
              <a class="text-emerald-400 hover:underline" href="cell-lines.html#${encodeURIComponent(it.cell_line_id)}">${escapeHtml(it.cell_line_id)}</a> ·
              <a class="text-emerald-400 hover:underline" href="modifiers.html#${encodeURIComponent(it.modifier_id)}">${escapeHtml(it.modifier_id)}</a>
            </div>

            <div class="mt-4 bg-slate-900/60 border border-slate-700 rounded-lg p-4 text-sm">
              <span class="kv-key">Evidence tier definition</span>
              <div class="text-slate-300 mt-1">${escapeHtml(tierDef.definition || "")}</div>
              ${it.combined_effect_metric !== null
                ? `<div class="mt-2 text-slate-300">Combined effect metric: <b class="text-white">${escapeHtml(it.combined_effect_metric)}</b></div>`
                : '<div class="mt-2 text-slate-500">No quantitative combined-effect metric recorded (qualitative claim).</div>'}
              ${it.model_run_id
                ? `<div class="mt-3 border border-violet-500/50 bg-violet-500/10 rounded-lg p-3 text-sm text-violet-200"><i class="fa-solid fa-flask mr-1"></i>Model-predicted interaction (Tier 2b) — model run <b>${escapeHtml(it.model_run_id)}</b>. Experimental: predictions are invalidated and regenerated as data and models improve, never silently overwritten.</div>`
                : ""}
            </div>

            <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-user-pen mr-2 text-emerald-400"></i>Curator notes</h3>
            <p class="text-slate-300 text-sm mt-1">${escapeHtml(it.curator_notes || "—")}</p>

            <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-dna mr-2 text-emerald-400"></i>Linked mechanism (stress signature)</h3>
            ${mechanismBlock(it.mechanism)}

            <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-book-open mr-2 text-emerald-400"></i>Sources</h3>
            <div class="mt-1">${citationChips(it.source_study)}</div>
            ${it.modifier?.protocol_parameters ? `
              <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-sliders mr-2 text-emerald-400"></i>Modifier protocol</h3>
              <div class="grid sm:grid-cols-2 md:grid-cols-3 gap-3 mt-2">${protocolKv(it.modifier.protocol_parameters)}</div>` : ""}
            ${relatedBlock(it.related || [])}
          </div>`;
    } catch (e) {
        box.innerHTML = `<div class="text-rose-400 text-sm">${escapeHtml(e.message)}</div>`;
    }
}

load();
