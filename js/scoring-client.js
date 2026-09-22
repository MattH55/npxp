/* Client-side synergy scoring — used by the static GitHub Pages build, where no
   backend exists. This is a line-faithful port of synlethality/scoring.py
   (Bliss / HSA / Loewe / ZIP over a modifier-level x drug-dose viability
   matrix); the only difference is that the Hill-curve least-squares fit uses a
   bounded Nelder-Mead simplex instead of scipy's curve_fit. Kept DOM-free so it
   can be unit-tested against the Python engine (see scripts/). */

const SCORING_MODELS = ["bliss", "hsa", "loewe", "zip"];
const SYNERGY_THRESHOLD = 0.1; // |mean excess| above this -> synergistic/antagonistic

/* ---------- input validation + effect matrix (E = 1 - viability) ---------- */

function asEffects(doses_a, doses_b, viability, viability_scale) {
    if (!Array.isArray(doses_a) || !Array.isArray(doses_b)) {
        throw new Error("doses_a and doses_b must be lists of numbers.");
    }
    if (!Array.isArray(viability) || !viability.length) {
        throw new Error("viability must be a non-empty 2-D list.");
    }
    if (!["fraction", "percent"].includes(viability_scale)) {
        throw new Error("viability_scale must be 'fraction' or 'percent'.");
    }
    const checkDoses = (d, name) => {
        if (d.length < 2) throw new Error(`${name} must contain a control (0) plus >=1 dose.`);
        const vals = d.map((x) => {
            if (typeof x !== "number" || Number.isNaN(x)) throw new Error(`${name} entries must be numeric, got ${JSON.stringify(x)}.`);
            return x;
        });
        if (vals[0] !== 0) throw new Error(`${name}[0] must be 0 (control edge of the matrix).`);
        for (let i = 1; i < vals.length; i++) {
            if (vals[i] <= vals[i - 1]) throw new Error(`${name} must be strictly ascending.`);
        }
        return vals;
    };
    const da = checkDoses(doses_a, "doses_a");
    const db = checkDoses(doses_b, "doses_b");
    const scale = viability_scale === "percent" ? 100 : 1;
    if (viability.length !== da.length) {
        throw new Error(`viability has ${viability.length} rows but doses_a has ${da.length} levels.`);
    }
    return [da, db, viability.map((row, i) => {
        if (!Array.isArray(row) || row.length !== db.length) {
            throw new Error(`viability row ${i} must be a list of ${db.length} numbers (one per doses_b entry).`);
        }
        return row.map((v) => {
            if (typeof v !== "number" || Number.isNaN(v)) throw new Error(`viability entries must be numeric, got ${JSON.stringify(v)}.`);
            const fv = v / scale;
            if (fv < -1e-9 || fv > 1 + 1e-9) throw new Error(`viability value ${JSON.stringify(v)} outside the ${viability_scale} range.`);
            return Math.min(1, Math.max(0, 1 - fv));
        });
    })];
}

/* ---------- Bliss / HSA (edge-null reference models) ---------- */

function edgeScores(effects, expectedFn) {
    return effects.map((row, i) => row.map((e, j) =>
        (i === 0 || j === 0) ? null : e - expectedFn(effects[i][0], effects[0][j])));
}
const blissScores = (effects) => edgeScores(effects, (ea, eb) => ea + eb - ea * eb);
const hsaScores = (effects) => edgeScores(effects, Math.max);

/* ---------- Hill-curve fitting (Loewe / ZIP) ---------- */

function hillEffect(d, emax, ec50, n) {
    if (d <= 0) return 0;
    return emax * Math.pow(d, n) / (Math.pow(ec50, n) + Math.pow(d, n));
}

function hillInverse(e, emax, ec50, n) {
    if (e <= 0 || e >= emax) return null;
    return ec50 * Math.pow(e / (emax - e), 1 / n);
}

function median(xs) {
    const s = [...xs].sort((a, b) => a - b);
    const m = Math.floor(s.length / 2);
    return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}


/* Bounded Nelder-Mead simplex (box constraints enforced by clamping). */
function nelderMead(f, x0, lb, ub, maxIter = 2000, tol = 1e-12) {
    const n = x0.length;
    const clamp = (x) => x.map((v, i) => Math.min(ub[i], Math.max(lb[i], v)));
    let simplex = [clamp(x0)];
    for (let i = 0; i < n; i++) {
        const x = [...simplex[0]];
        x[i] += x[i] !== 0 ? x[i] * 0.05 : (ub[i] - lb[i]) * 0.01;
        simplex.push(clamp(x));
    }
    let vals = simplex.map(f);
    const sort = () => {
        const idx = simplex.map((_, i) => i).sort((a, b) => vals[a] - vals[b]);
        simplex = idx.map((i) => simplex[i]);
        vals = idx.map((i) => vals[i]);
    };
    for (let iter = 0; iter < maxIter; iter++) {
        sort();
        if (Math.abs(vals[n] - vals[0]) <= tol * (Math.abs(vals[0]) + tol)) break;
        const centroid = Array(n).fill(0);
        for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) centroid[k] += simplex[i][k] / n;
        const worst = simplex[n];
        const refl = clamp(centroid.map((c, k) => c + (c - worst[k])));
        const fR = f(refl);
        if (fR < vals[0]) {
            const exp = clamp(centroid.map((c, k) => c + 2 * (c - worst[k])));
            const fE = f(exp);
            if (fE < fR) { simplex[n] = exp; vals[n] = fE; } else { simplex[n] = refl; vals[n] = fR; }
        } else if (fR < vals[n - 1]) {
            simplex[n] = refl; vals[n] = fR;
        } else {
            const con = clamp(centroid.map((c, k) => fR < vals[n] ? c + 0.5 * (refl[k] - c) : c + 0.5 * (worst[k] - c)));
            const fC = f(con);
            if (fC < Math.min(fR, vals[n])) {
                simplex[n] = con; vals[n] = fC;
            } else {
                for (let i = 1; i <= n; i++) {
                    simplex[i] = clamp(simplex[0].map((b, k) => b + 0.5 * (simplex[i][k] - b)));
                    vals[i] = f(simplex[i]);
                }
            }
        }
    }
    sort();
    return simplex[0];
}

/* Fit E(d) = Emax * d^n / (EC50^n + d^n), Emin fixed at 0 — same p0/bounds and
   the same identifiability rule as the Python engine: with only 3 non-control
   points the Hill slope is fixed at n = 1; >=4 points fit the full curve. */
function fitHill(doses, edgeEffects) {
    const pts = doses.map((d, i) => [d, edgeEffects[i]]).filter(([d]) => d > 0);
    if (pts.length < 3) {
        throw new Error(`Loewe/ZIP need >=3 non-control dose points on each monotherapy edge; got ${pts.length}.`);
    }
    const d = pts.map((p) => p[0]), e = pts.map((p) => p[1]);
    const emax0 = Math.min(1, Math.max(Math.max(...e), 0.05));
    const dMin = Math.min(...d), dMax = Math.max(...d), med = median(d);
    const fitSlope = pts.length >= 4;
    const p0 = fitSlope ? [emax0, med, 1.0] : [emax0, med];
    const lb = fitSlope ? [1e-4, dMin * 1e-2, 0.1] : [1e-4, dMin * 1e-2];
    const ub = fitSlope ? [1.0, dMax * 1e2, 10.0] : [1.0, dMax * 1e2];
    const sse = (p) => {
        const n = fitSlope ? p[2] : 1.0;
        let s = 0;
        for (let i = 0; i < d.length; i++) {
            const r = hillEffect(d[i], p[0], p[1], n) - e[i];
            s += r * r;
        }
        return s;
    };
    const p = nelderMead(sse, p0, lb, ub);
    return [p[0], p[1], fitSlope ? p[2] : 1.0];
}

/* ---------- Loewe / ZIP ---------- */

function loeweScores(doses_a, doses_b, effects) {
    const hillA = fitHill(doses_a, effects.map((row) => row[0]));
    const hillB = fitHill(doses_b, effects[0]);
    const cis = [];
    const out = effects.map((row, i) => row.map((e, j) => {
        if (i === 0 || j === 0) return null;
        const dA = hillInverse(e, ...hillA);
        const dB = hillInverse(e, ...hillB);
        if (dA === null || dB === null) return null; // effect outside what either agent can explain
        const ci = doses_a[i] / dA + doses_b[j] / dB;
        cis.push(ci);
        return 1 - ci;
    }));
    return [out, {
        hill_a: hillA, hill_b: hillB,
        mean_ci: cis.length ? cis.reduce((a, b) => a + b, 0) / cis.length : null,
    }];
}

function zipScores(doses_a, doses_b, effects) {
    const hillA = fitHill(doses_a, effects.map((row) => row[0]));
    const hillB = fitHill(doses_b, effects[0]);
    const out = effects.map((row, i) => row.map((e, j) => {
        if (i === 0 || j === 0) return null;
        const fa = hillEffect(doses_a[i], ...hillA);
        const fb = hillEffect(doses_b[j], ...hillB);
        return e - (fa + fb - fa * fb);
    }));
    return [out, { hill_a: hillA, hill_b: hillB }];
}

/* ---------- dispatcher ---------- */

function summarize(scores) {
    const vals = scores.flat().filter((v) => v !== null);
    if (!vals.length) return { n_cells: 0, mean: null, max: null, min: null, classification: "no_data" };
    const mean = vals.reduce((a, b) => a + b, 0) / vals.length;
    const classification = mean > SYNERGY_THRESHOLD ? "synergistic"
        : mean < -SYNERGY_THRESHOLD ? "antagonistic" : "additive";
    return { n_cells: vals.length, mean, max: Math.max(...vals), min: Math.min(...vals), classification };
}

function scoreMatrixClient(doses_a, doses_b, viability, models = null, viability_scale = "fraction") {
    const wanted = models && models.length ? models : SCORING_MODELS;
    const unknown = wanted.filter((m) => !SCORING_MODELS.includes(m));
    if (unknown.length) throw new Error(`Unknown scoring model(s) ${JSON.stringify(unknown)}; choose from ${JSON.stringify(SCORING_MODELS)}.`);
    const [da, db, effects] = asEffects(doses_a, doses_b, viability, viability_scale);
    const results = {};
    for (const name of wanted) {
        let scores, meta = null;
        if (name === "bliss") scores = blissScores(effects);
        else if (name === "hsa") scores = hsaScores(effects);
        else if (name === "loewe") [scores, meta] = loeweScores(da, db, effects);
        else [scores, meta] = zipScores(da, db, effects);
        results[name] = { scores, summary: summarize(scores), meta };
    }
    return {
        models: results,
        n_modifier_levels: da.length,
        n_drug_doses: db.length,
        synergy_threshold: SYNERGY_THRESHOLD,
    };
}

/* Node export for parity tests against the Python engine; ignored by browsers. */
if (typeof module !== "undefined" && module.exports) {
    module.exports = { scoreMatrixClient, fitHill, nelderMead, asEffects, SCORING_MODELS, SYNERGY_THRESHOLD };
}
