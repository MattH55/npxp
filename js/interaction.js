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

function clinicalEvidenceBlock(rows) {
    if (!rows || !rows.length) {
        return `<div class="text-slate-500 text-sm">No registered human trial links this specific
          modifier × drug × cell-line finding yet — that is not evidence the combination has
          never been tested, only that no match has been curated here.</div>`;
    }
    return rows.map((c) => `
      <div class="bg-slate-900/60 border border-slate-700 rounded-lg p-4 text-sm mb-3">
        <div class="flex flex-wrap items-center gap-2">
          <a class="cite-chip" target="_blank" rel="noopener" href="https://clinicaltrials.gov/study/${encodeURIComponent(c.nct_id)}">${escapeHtml(c.nct_id)} <i class="fa-solid fa-arrow-up-right-from-square text-[0.6rem]"></i></a>
          ${c.phase ? `<span class="cite-chip">${escapeHtml(c.phase)}</span>` : ""}
          ${c.status ? `<span class="cite-chip">${escapeHtml(c.status)}</span>` : ""}
        </div>
        <p class="text-slate-300 mt-2">${escapeHtml(c.outcome_summary || "")}</p>
        ${c.publication_ref ? `<div class="mt-2">${citationChips([c.publication_ref])}</div>` : ""}
      </div>`).join("");
}

function relatedBlock(related, truncated, cellLineId) {
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
            <td class="py-2 pr-3 text-white">${escapeHtml(r.drug_name)}${r.drug_lincs_signature_ref ? ' <i class="fa-solid fa-dna text-[0.6rem] text-emerald-400" title="LINCS-backed"></i>' : ""}</td>
            <td class="py-2 pr-3 text-slate-300">${escapeHtml(r.modifier_agent)}</td>
            <td class="py-2 pr-3">${typeBadge(r.interaction_type)}</td>
            <td class="py-2 pr-3">${tierBadge(r.evidence_tier)}</td>
          </tr>`).join("")}
        </tbody></table>
        ${truncated ? `<div class="text-slate-500 text-xs mt-2">Showing the first ${related.length} (Tier 1-3 only) — this cell line has more. See <a class="underline text-emerald-400" href="cell-lines.html#${encodeURIComponent(cellLineId)}">the cell line page</a> or filter <code>/api/interactions?cell_line_id=${encodeURIComponent(cellLineId)}</code> for the rest, including Tier 4.</div>` : ""}`;
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
            <div class="flex flex-wrap items-center justify-between gap-2 mb-3">
              <div class="flex flex-wrap items-center gap-2">
                ${typeBadge(it.interaction_type)} ${tierBadge(it.evidence_tier)}
                <span class="text-slate-500 text-xs ml-1">${escapeHtml(tierDef.label || "")}</span>
                ${lowResourceBadge(it.low_resource_relevance)}
              </div>
              <button disabled title="Pending confirmation of the trial builder tool's integration interface (build spec Open Questions) — not wired up yet."
                      class="bg-slate-800 border border-slate-700 text-slate-500 text-xs font-semibold px-3 py-1.5 rounded-lg cursor-not-allowed">
                <i class="fa-solid fa-paper-plane mr-1"></i>Send to Trial Builder
              </button>
            </div>
            <h1 class="text-2xl font-extrabold text-white leading-snug">
              <a class="hover:underline" href="drugs.html#${encodeURIComponent(it.drug_id)}">${escapeHtml(it.drug_name)}</a>
              <span class="text-slate-500 font-normal">+</span> ${escapeHtml(it.modifier_agent)}
              <span class="text-slate-500 font-normal">in</span> ${escapeHtml(it.cell_line_name)}
            </h1>
            <div class="text-slate-400 text-sm mt-1">${escapeHtml(it.tissue_origin || "")} ·
              <a class="text-emerald-400 hover:underline" href="cell-lines.html#${encodeURIComponent(it.cell_line_id)}">${escapeHtml(it.cell_line_id)}</a> ·
              <a class="text-emerald-400 hover:underline" href="modifiers.html#${encodeURIComponent(it.modifier_id)}">${escapeHtml(it.modifier_id)}</a> ·
              <a class="text-emerald-400 hover:underline" href="drugs.html#${encodeURIComponent(it.drug_id)}">${escapeHtml(it.drug_id)}</a>
            </div>
            <div class="text-xs mt-2 flex flex-wrap gap-2 items-center">
              ${it.drug_lincs_signature_ref
                ? `<span class="cite-chip">LINCS signature: ${escapeHtml(it.drug_lincs_signature_ref)}</span>`
                : '<span class="text-slate-500">No LINCS L1000 signature for this drug yet — not eligible as Step 2 model input until ingested.</span>'}
              ${it.drug_pubchem_cid ? `<a class="cite-chip" target="_blank" rel="noopener" href="https://pubchem.ncbi.nlm.nih.gov/compound/${encodeURIComponent(it.drug_pubchem_cid)}">PubChem ${escapeHtml(it.drug_pubchem_cid)} <i class="fa-solid fa-arrow-up-right-from-square text-[0.6rem]"></i></a>` : ""}
            </div>

            <div class="mt-4 bg-slate-900/60 border border-slate-700 rounded-lg p-4 text-sm">
              <span class="kv-key">Evidence tier definition</span>
              <div class="text-slate-300 mt-1">${escapeHtml(tierDef.definition || "")}</div>
              ${it.combined_effect_metric !== null
                ? `<div class="mt-2 text-slate-300">Combined effect metric: <b class="text-white">${escapeHtml(it.combined_effect_metric)}</b></div>`
                : '<div class="mt-2 text-slate-500">No quantitative combined-effect metric recorded (qualitative claim).</div>'}
              ${it.synergy_model_scores
                ? `<div class="mt-2 text-slate-300">
                     <div class="text-xs text-slate-500 mb-1">Per-model synergy scores (Bliss/HSA/Loewe/ZIP) — shown separately, not collapsed to one number, since these can genuinely disagree (see <a class="underline" href="candidates.html">methodology</a>):</div>
                     <div class="flex flex-wrap gap-2">
                       ${Object.entries(it.synergy_model_scores).map(([name, d]) => `
                         <span class="cite-chip" title="${d.n_cells} cell(s) scored">
                           ${name.toUpperCase()}: <b>${d.mean !== null ? d.mean.toFixed(3) : "n/a"}</b> (${escapeHtml(d.classification)})
                         </span>`).join("")}
                     </div>
                   </div>`
                : ""}
              <div class="mt-2 text-slate-400">Priority score: <b class="text-white">${it.priority_score}</b>
                <a class="text-emerald-400 hover:underline ml-1" href="candidates.html">(methodology)</a></div>
              ${it.model_run_id
                ? `<div class="mt-3 border border-violet-500/50 bg-violet-500/10 rounded-lg p-3 text-sm text-violet-200"><i class="fa-solid fa-flask mr-1"></i>Model-predicted interaction (Tier 2b) — model run <b>${escapeHtml(it.model_run_id)}</b>. Experimental: predictions are invalidated and regenerated as data and models improve, never silently overwritten.</div>`
                : ""}
              ${it.is_bridging_candidate
                ? `<div class="mt-3 border border-violet-500/50 bg-violet-500/10 rounded-lg p-3 text-sm text-violet-200"><i class="fa-solid fa-flask-vial mr-1"></i>Mechanistic Bridging output — a <b>model-generated hypothesis</b> from matching this modifier's induced signature against ${escapeHtml(it.drug_name)}'s known resistance mechanism, not a tested result. See <a class="underline" href="candidates.html">Candidates</a> for the full list.</div>`
                : ""}
              ${it.evidence_tier === "tier_4_heuristic_target_match"
                ? `<div class="mt-3 border border-red-500/50 bg-red-500/10 rounded-lg p-3 text-sm text-red-200"><i class="fa-solid fa-triangle-exclamation mr-1"></i><b>Tier 4 — mechanical keyword match, not a verified claim.</b> ${escapeHtml(it.drug_name)}'s own DepMap-supplied target/mechanism text happens to share a keyword with this modifier's induced signature category. This is a purely textual co-occurrence, not a cited resistance mechanism (compare Tier 3), has not been reviewed by a curator, and predicts no direction. Treat as a keyword-matched lead worth checking, nothing more.</div>`
                : ""}
            </div>

            <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-user-pen mr-2 text-emerald-400"></i>Curator notes</h3>
            <p class="text-slate-300 text-sm mt-1">${escapeHtml(it.curator_notes || "—")}</p>

            <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-dna mr-2 text-emerald-400"></i>Linked mechanism (stress signature)</h3>
            ${mechanismBlock(it.mechanism)}
            ${it.resistance_mechanism ? `
              <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-shield-halved mr-2 text-emerald-400"></i>Drug resistance mechanism matched</h3>
              <table class="w-full text-sm mt-2"><tbody>
                <tr class="border-b border-slate-700/60"><td class="py-2 pr-3 kv-key w-40">Pathway/gene</td><td class="py-2 text-slate-200">${escapeHtml(it.resistance_mechanism.pathway_or_gene)}</td></tr>
                <tr><td class="py-2 pr-3 kv-key">Description</td><td class="py-2 text-slate-300">${escapeHtml(it.resistance_mechanism.mechanism_description)}</td></tr>
              </tbody></table>` : ""}

            <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-hospital mr-2 text-emerald-400"></i>Human clinical evidence</h3>
            ${clinicalEvidenceBlock(it.clinical_evidence)}

            <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-book-open mr-2 text-emerald-400"></i>Sources</h3>
            <div class="mt-1">${citationChips(it.source_study)}</div>
            ${it.modifier?.protocol_parameters ? `
              <h3 class="text-white font-semibold mt-6"><i class="fa-solid fa-sliders mr-2 text-emerald-400"></i>Modifier protocol</h3>
              <div class="grid sm:grid-cols-2 md:grid-cols-3 gap-3 mt-2">${protocolKv(it.modifier.protocol_parameters)}</div>` : ""}
            ${relatedBlock(it.related || [], it.related_truncated, it.cell_line_id)}
          </div>`;
    } catch (e) {
        box.innerHTML = `<div class="text-rose-400 text-sm">${escapeHtml(e.message)}</div>`;
    }
}

load();
