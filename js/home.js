/* Home page: stats + tier legend + curated sources. */
renderNav("Home");

const SOURCES = [
    "Zhang Z et al., Oncol Lett 2021 — glucose deprivation & β-hydroxybutyrate RNA-seq, MCF-7/T47D (GSE153830; PMID 33281976)",
    "GEO GSE48398 — hyperthermic shock (42–45 °C), breast cancer panel vs MCF-10A (PMID 27245201)",
    "Tabuchi Y et al., Int J Hyperthermia 2008 — mild hyperthermia 41 °C/30 min, U-937; HSF1/Hsp40/Hsp70, no apoptosis (GSE10043; PMID 18608577)",
    "Helderman RFCPA et al., Cells 2020 — HIPEC-mimetic hyperthermia × platinum drugs/5-FU/MMC in CRC lines (PMID 32722384)",
    "Raffaele M et al., Int J Mol Sci 2019 — metformin × glucose deprivation, DU-145; HO-1 (PMID 31137785)",
    "Rohwer N & Cramer T, Drug Resist Updat 2011 — HIF-1 hypoxia resistance review (PMID 21466972)",
    "Hadad SM et al., Mol Med Rep 2014 — metformin impairs breast cancer cell growth (PMID 24338509)",
    "Yang J et al., Int J Biol Sci 2021 — metformin induces ferroptosis via SLC7A11 (PMID 34162423)",
];
document.getElementById("source-list").innerHTML =
    SOURCES.map((s) => `<li>${escapeHtml(s)}</li>`).join("");

(async () => {
    try {
        const [stats, tiers] = await Promise.all([apiGet("/stats"), apiGet("/evidence_tiers")]);
        document.getElementById("stat-cell-lines").textContent = stats.cell_lines;
        document.getElementById("stat-modifiers").textContent = stats.modifiers;
        document.getElementById("stat-drugs").textContent = stats.drugs;
        document.getElementById("stat-interactions").textContent = stats.interaction_effects;

        const tierOrder = ["tier_1_direct", "tier_2_inferred", "tier_2b_model_predicted", "tier_3_mechanism_only"];
        const tierColors = {
            tier_1_direct: "bg-emerald-500",
            tier_2_inferred: "bg-yellow-500",
            tier_2b_model_predicted: "bg-violet-500",
            tier_3_mechanism_only: "bg-orange-500",
        };
        const maxTier = Math.max(1, ...Object.values(stats.interactions_by_tier));
        document.getElementById("tier-bars").innerHTML = tierOrder.map((t) => {
            const n = stats.interactions_by_tier[t] || 0;
            const pct = Math.round((n / maxTier) * 100);
            return `<div>
              <div class="flex justify-between mb-1"><span>${escapeHtml(tiers[t]?.label || t)}</span><span class="text-slate-400">${n}</span></div>
              <div class="w-full bg-slate-700 rounded h-2"><div class="${tierColors[t]} h-2 rounded" style="width:${pct}%"></div></div>
            </div>`;
        }).join("");

        const typeOrder = ["synergistic", "additive", "antagonistic", "unknown"];
        const maxType = Math.max(1, ...Object.values(stats.interactions_by_type));
        document.getElementById("type-bars").innerHTML = typeOrder.map((t) => {
            const n = stats.interactions_by_type[t] || 0;
            const pct = Math.round((n / maxType) * 100);
            return `<div>
              <div class="flex justify-between mb-1">${typeBadge(t)}<span class="text-slate-400">${n}</span></div>
              <div class="w-full bg-slate-700 rounded h-2"><div class="bg-slate-400 h-2 rounded" style="width:${pct}%"></div></div>
            </div>`;
        }).join("");

        document.getElementById("tier-cards").innerHTML = tierOrder.map((t) => {
            const d = tiers[t];
            return `<div class="bg-slate-800/60 border border-slate-700 rounded-xl p-5">
              <div class="mb-2">${tierBadge(t)} <span class="text-white font-semibold ml-1">${escapeHtml(d.label)}</span></div>
              <p class="text-slate-400 text-sm">${escapeHtml(d.definition)}</p>
            </div>`;
        }).join("");
    } catch (e) {
        document.getElementById("stat-cards").insertAdjacentHTML(
            "afterend", `<div class="text-rose-400 text-sm mb-6">Failed to load stats: ${escapeHtml(e.message)}</div>`);
    }
})();
