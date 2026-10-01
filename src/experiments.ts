// Experiment suite. Writes results/results.json (consumed by scripts/plot.py)
// and results/results.md. Use --quick for a fast smoke run.

import { mkdirSync, writeFileSync } from 'node:fs';
import { performance } from 'node:perf_hooks';
import { CableModel, PASSIVE, Simulator, type Biophysics, type Method } from './cable.ts';
import { countFill, denseSolve, orders } from './dense.ts';
import { buildTree, cable, rallTree, randomNodeTree, randomSections, rng, type BranchMode, type SectionSpec } from './morphology.ts';
import { equivalentCylinder, lambdaCm, rmFromLeak, treeInputResistance } from './analytic.ts';
import { HH_LEAK, HH_SQUID } from './hh.ts';

const QUICK = process.argv.includes('--quick');
const out: Record<string, unknown> = {};
const md: string[] = ['# Results', '', QUICK ? '_Quick mode: reduced sizes._' : '_Full run._', ''];
const log = (s = '') => {
  console.log(s);
  md.push(s);
};
const fmt = (x: number, p = 3) =>
  Number.isInteger(x) && Math.abs(x) < 1e7
    ? String(x)
    : Math.abs(x) >= 1e-3 && Math.abs(x) < 1e5
      ? x.toPrecision(p)
      : x.toExponential(p - 1);
const table = (head: string[], rows: (string | number)[][]) => {
  log('| ' + head.join(' | ') + ' |');
  log('|' + head.map(() => '---').join('|') + '|');
  for (const r of rows) log('| ' + r.map((c) => (typeof c === 'number' ? fmt(c) : c)).join(' | ') + ' |');
  log();
};
/** Least-squares slope of log(y) against log(x). */
const slope = (xs: number[], ys: number[]) => {
  const lx = xs.map(Math.log), ly = ys.map(Math.log);
  const mx = lx.reduce((a, b) => a + b) / lx.length, my = ly.reduce((a, b) => a + b) / ly.length;
  let num = 0, den = 0;
  lx.forEach((x, i) => {
    num += (x - mx) * (ly[i]! - my);
    den += (x - mx) ** 2;
  });
  return num / den;
};
const median = (xs: number[]) => [...xs].sort((a, b) => a - b)[Math.floor(xs.length / 2)]!;

const RM = rmFromLeak(PASSIVE.gLeak);
const tree = (secs: SectionSpec[], branch: BranchMode) => buildTree(secs, { branch, rootNode: true });

// ---------------------------------------------------------------------------
// E1. The solver: correctness, fill-in, scaling.
// ---------------------------------------------------------------------------
function e1() {
  log('## E1. Hines solver: correctness, fill-in, and scaling');
  log();
  // Correctness on real backward-Euler cable matrices from random dendrites.
  let worst = 0;
  for (let seed = 1; seed <= (QUICK ? 5 : 30); seed++) {
    const t = buildTree(randomSections(8 + (seed % 12), 1 + (seed % 9), seed), { branch: 'junction', rootNode: true });
    const model = new CableModel(t.morph, PASSIVE);
    const A = model.matrix(model.cap.map((c, i) => c / 0.025 + model.area[i]! * PASSIVE.gLeak));
    const r = rng(seed);
    const b = Float64Array.from({ length: model.n }, () => r() - 0.5);
    const ref = denseSolve(A.toDense(), b);
    const x = A.solve(Float64Array.from(b));
    const scale = Math.max(...ref.map(Math.abs));
    for (let i = 0; i < model.n; i++) worst = Math.max(worst, Math.abs(x[i]! - ref[i]!) / scale);
  }
  log(`Max relative difference from dense LU with partial pivoting on ${QUICK ? 5 : 30} backward-Euler cable matrices (random dendrites with junction nodes): **${fmt(worst)}**.`);
  log();

  // Fill-in under different elimination orders.
  const fillRows: (string | number)[][] = [];
  const fill: Record<string, number[]> = { n: [], hines: [], rootFirst: [], random: [] };
  for (const n of QUICK ? [100, 300] : [100, 300, 1000]) {
    const m = randomNodeTree(n, n);
    const o = orders(n, rng(n + 1));
    const f = { hines: countFill(m.parent, o.hines), rootFirst: countFill(m.parent, o.rootFirst), random: countFill(m.parent, o.random) };
    fill.n!.push(n);
    fill.hines!.push(f.hines);
    fill.rootFirst!.push(f.rootFirst);
    fill.random!.push(f.random);
    fillRows.push([n, n - 1, f.hines, f.rootFirst, f.random, (n * (n - 1)) / 2]);
  }
  log('Fill-in (new off-diagonal pairs) from symbolic elimination on uniform random recursive trees:');
  log();
  table(['n', 'tree edges', 'leaves-first (Hines)', 'root-first', 'random order', 'dense n(n-1)/2'], fillRows);

  // Wall-clock scaling.
  const hinesSizes = QUICK ? [1e3, 1e4, 1e5] : [1e3, 1e4, 1e5, 1e6];
  const denseSizes = QUICK ? [100, 200] : [100, 200, 400, 800];
  const time = (f: () => void, reps: number) => {
    f();
    const ts: number[] = [];
    for (let k = 0; k < reps; k++) {
      const t0 = performance.now();
      f();
      ts.push(performance.now() - t0);
    }
    return median(ts);
  };
  const build = (n: number) => {
    const m = randomNodeTree(n, 99);
    const model = new CableModel(m, PASSIVE);
    return model.matrix(model.cap.map((c) => c / 0.025));
  };
  const hinesMs = hinesSizes.map((n) => {
    const A = build(n);
    const b = new Float64Array(n).fill(1);
    const x = new Float64Array(n);
    return time(() => {
      x.set(b);
      A.solve(x);
    }, n >= 1e6 ? 15 : 51);
  });
  const denseMs = denseSizes.map((n) => {
    const A = build(n).toDense();
    const b = new Float64Array(n).fill(1);
    return time(() => denseSolve(A, b), n >= 800 ? 3 : 7);
  });
  const hinesAt = (n: number) => {
    const A = build(n);
    const x = new Float64Array(n);
    return time(() => {
      x.fill(1);
      A.solve(x);
    }, 101);
  };
  log('Wall-clock time per solve (Node ' + process.version + ', single thread, median):');
  log();
  table(['n', 'Hines (ms)', 'ns / node'], hinesSizes.map((n, i) => [n, hinesMs[i]!, (hinesMs[i]! * 1e6) / n]));
  table(['n', 'dense LU (ms)', 'Hines at same n (ms)', 'speedup'], denseSizes.map((n, i) => {
    const h = hinesAt(n);
    return [n, denseMs[i]!, h, denseMs[i]! / h];
  }));
  log(`Fitted exponents: Hines **${fmt(slope(hinesSizes, hinesMs))}**, dense LU **${fmt(slope(denseSizes, denseMs))}**.`);
  log();
  out.e1 = { worst, fill, hinesSizes, hinesMs, denseSizes, denseMs };
}

// ---------------------------------------------------------------------------
// E2. Spatial convergence against exact cable theory; Rall equivalence.
// ---------------------------------------------------------------------------
function e2() {
  log('## E2. Spatial convergence against exact cable theory');
  log();
  const nsegs = QUICK ? [1, 2, 4, 8, 16] : [1, 2, 4, 8, 16, 32, 64, 128];
  const geoms: Record<string, (nseg: number) => SectionSpec[]> = {
    cable: (nseg) => [{ parent: -1, length: 1000, diam: 2, nseg }],
    rall: (nseg) => rallTree(4, 300, 4, nseg),
    random: (nseg) => randomSections(25, nseg, 5),
  };
  const conv: Record<string, Record<string, number[]>> = {};
  const rows: (string | number)[][] = [];
  for (const [name, g] of Object.entries(geoms)) {
    const exact = treeInputResistance(g(1), RM, PASSIVE.ra);
    for (const mode of ['junction', 'direct'] as const) {
      if (name === 'cable' && mode === 'direct') continue; // identical without branch points
      const errs = nsegs.map((k) => {
        const r = new CableModel(tree(g(k), mode).morph, PASSIVE).inputResistance(0);
        return Math.abs(r - exact) / exact;
      });
      conv[`${name}/${mode}`] = { nsegs, errs };
      const tail = 4;
      rows.push([name, mode, fmt(exact, 6), errs[0]!, errs[errs.length - 1]!, -slope(nsegs.slice(-tail), errs.slice(-tail))]);
    }
  }
  log('Relative error of the steady-state input resistance at the root, compared with the exact');
  log('leaves-up recursion (Rm = 10 kΩ·cm², Ra = 100 Ω·cm). The order is fitted over the four finest grids.');
  log();
  table(['geometry', 'branch points', 'exact R_in (MΩ)', `err @ nseg=${nsegs[0]}`, `err @ nseg=${nsegs.at(-1)}`, 'observed order'], rows);

  // Rall equivalence, steady state and transient.
  const rall = rallTree(4, 300, 4, 1);
  const eq = equivalentCylinder(rall, RM, PASSIVE.ra);
  const rTree = treeInputResistance(rall, RM, PASSIVE.ra);
  const rEq = treeInputResistance([eq], RM, PASSIVE.ra);
  const broken = rall.map((s, i) => (i === 0 ? s : { ...s, diam: s.diam * 1.3 }));
  // Backward Euler: a current step excites every spatial mode, and CN would
  // ring on the stiff ones (see E3). The same dt is used for tree and cylinder.
  const dt = 0.002;
  const T = QUICK ? 10 : 30;
  const every = Math.round(0.05 / dt);
  const runStep = (secs: SectionSpec[], nseg: number) => {
    // Same number of compartments per λ in every section.
    const lam0 = lambdaCm(secs[0]!.diam, RM, PASSIVE.ra);
    const t = tree(secs.map((s) => {
      const X = (s.length * 1e-4) / lambdaCm(s.diam, RM, PASSIVE.ra);
      return { ...s, nseg: Math.max(1, Math.round((nseg * X * lam0) / (secs[0]!.length * 1e-4))) };
    }), 'junction');
    const sim = new Simulator(new CableModel(t.morph, PASSIVE), 'be', dt);
    sim.addStim({ node: 0, amp: 0.1, t0: 0, t1: 1e9 });
    const v: number[] = [0];
    let k = 0;
    sim.run(T, (s) => {
      if (++k % every === 0) v.push(s.V[0]!);
    });
    return v;
  };
  const maxDiff = (a: number[], b: number[]) => Math.max(...a.map((v, i) => Math.abs(v - b[i]!))) / Math.max(...b);
  const rallRows: (string | number)[][] = [];
  const rallConv: { nseg: number; dTree: number; dBroken: number }[] = [];
  let traces: Record<string, number[]> = {};
  for (const nseg of QUICK ? [4, 8] : [4, 8, 16, 32]) {
    const cyl = runStep([{ ...eq, nseg: 1 }], nseg);
    const tr = runStep(rall, nseg);
    const br = runStep(broken, nseg);
    const dTree = maxDiff(tr, cyl), dBroken = maxDiff(br, cyl);
    rallConv.push({ nseg, dTree, dBroken });
    rallRows.push([nseg, fmt(dTree * 100) + '%', fmt(dBroken * 100) + '%']);
    traces = { tree: tr, cylinder: cyl, broken: br };
  }
  log(`**Rall's equivalent cylinder.** A 4-generation binary tree obeying d_p^{3/2} = Σ d_c^{3/2}, with equal electrotonic path lengths, has exact R_in = ${fmt(rTree, 8)} MΩ. Its equivalent cylinder (L = ${fmt(eq.length, 5)} µm) has ${fmt(rEq, 8)} MΩ.`);
  log(`The transient check: under a 0.1 nA step, compare the simulated root voltage of the tree and the cylinder (BE, dt = ${dt} ms, ${T} ms), as max |ΔV| / plateau:`);
  log();
  table(['compartments in root section', 'tree vs cylinder', '3/2 rule broken (daughters ×1.3 diameter)'], rallRows);
  log('With the 3/2 rule, the tree–cylinder gap is a discretisation artifact that shrinks under refinement. When the rule is broken, a real, grid-independent difference remains.');
  log();
  out.e2 = { conv, rall: { rTree, rEq, rallConv, traces, sampleDt: dt * every } };
}

// ---------------------------------------------------------------------------
// E3. Temporal convergence (BE vs staggered CN) and Crank–Nicolson ringing.
// ---------------------------------------------------------------------------
const HH_BIO: Biophysics = { cm: 1, ra: 100, ...HH_LEAK, hh: HH_SQUID };

function e3() {
  log('## E3. Temporal convergence: backward Euler vs staggered Crank–Nicolson');
  log();
  const secs = rallTree(3, 200, 4, 10);
  const t = tree(secs, 'junction');
  const probe = t.last[secs.length - 1]!;
  const T = QUICK ? 8 : 15;
  const sample = 0.1;
  const runTrace = (bio: Biophysics, method: Method, dt: number, amp: number) => {
    const sim = new Simulator(new CableModel(t.morph, bio), method, dt, bio.hh ? -65 : bio.eLeak);
    sim.addStim({ node: 0, amp, t0: 1, t1: 3 });
    const every = Math.round(sample / dt);
    const v: number[] = [sim.V[probe]!];
    let k = 0;
    sim.run(T, (s) => {
      if (++k % every === 0) v.push(s.V[probe]!);
    });
    return v;
  };
  const dts = QUICK ? [0.1, 0.05, 0.025] : [0.1, 0.05, 0.025, 0.0125, 0.00625];
  const refDt = QUICK ? 0.1 / 64 : 0.1 / 512;
  const res: Record<string, Record<string, number[]>> = {};
  const rows: (string | number)[][] = [];
  const hhTraces: Record<string, number[]> = {};
  for (const [label, bio, amp] of [['passive', PASSIVE, 0.2], ['HH', HH_BIO, 0.3]] as const) {
    const ref = runTrace(bio, 'cn', refDt, amp);
    if (label === 'HH') hhTraces.ref = ref;
    res[label] = { dts };
    for (const method of ['be', 'cn'] as const) {
      const errs = dts.map((dt) => {
        const v = runTrace(bio, method, dt, amp);
        if (label === 'HH' && dt === 0.1) hhTraces[method] = v;
        return Math.max(...v.map((x, i) => Math.abs(x - ref[i]!)));
      });
      res[label]![method] = errs;
      rows.push([label, method.toUpperCase(), errs[0]!, errs.at(-1)!, slope(dts, errs)]);
    }
  }
  log(`Max |V − V_ref| (mV) at a distal tip of a 3-generation tree, sampled every ${sample} ms over ${T} ms (2 ms current step).`);
  log(`Reference: CN at dt = ${fmt(refDt)} ms. HH is the squid axon at 6.3 °C on every compartment.`);
  log();
  table(['membrane', 'method', `err @ dt=${dts[0]}`, `err @ dt=${dts.at(-1)}`, 'observed order'], rows);

  // Ringing: CN is A-stable but its amplification factor tends to -1 for stiff modes.
  const nseg = 200;
  const ringT = 10;
  const ringModel = new CableModel(cable(1000, 2, nseg, true).morph, PASSIVE);
  const lam = new CableModel(cable(1000, 2, nseg, false).morph, PASSIVE).maxRate();
  const stepTrace = (method: Method, dt: number, sampleEvery: number) => {
    const sim = new Simulator(ringModel, method, dt);
    sim.addStim({ node: 0, amp: 0.1, t0: 0, t1: 1e9 });
    const v: number[] = [0];
    let k = 0;
    sim.run(ringT, (s) => {
      if (++k % sampleEvery === 0) v.push(s.V[0]!);
    });
    return v;
  };
  // Reference: CN with λ_max·dt < 2, so every mode decays monotonically. Sampled every 0.05 ms.
  const refStep = 0.05 / Math.ceil(0.05 / (1 / lam));
  const refTrace = stepTrace('cn', refStep, Math.round(0.05 / refStep));
  const plateau = ringModel.inputResistance(0) * 0.1;
  const ring: Record<string, { dt: number; v: number[] }> = { ref: { dt: 0.05, v: refTrace } };
  const ringRows: (string | number)[][] = [];
  for (const dt of [0.5, 0.05]) {
    for (const method of ['be', 'cn'] as const) {
      const v = stepTrace(method, dt, 1);
      const stride = Math.round(dt / 0.05);
      let reversals = 0;
      for (let i = 2; i < v.length; i++) if ((v[i]! - v[i - 1]!) * (v[i - 1]! - v[i - 2]!) < 0) reversals++;
      const err = v.map((x, i) => Math.abs(x - refTrace[i * stride]!) / plateau);
      ring[`${method}@${dt}`] = { dt, v };
      ringRows.push([method.toUpperCase(), dt, reversals, fmt(Math.max(...err.slice(1)) * 100) + '%', fmt(err.at(-1)! * 100) + '%']);
    }
  }
  log(`**Crank–Nicolson ringing.** Root voltage of a sealed 1000 µm cable (nseg = ${nseg}) under a 0.1 nA step. The exact response rises monotonically.`);
  log(`The stiffest mode has λ_max ≈ ${fmt(lam)} /ms. At dt = 0.5 ms its CN amplification factor is (1 − λdt/2)/(1 + λdt/2) ≈ ${fmt((1 - (lam * 0.5) / 2) / (1 + (lam * 0.5) / 2), 4)}, so it flips sign every step and barely decays.`);
  log(`The reference is CN at dt = ${fmt(refStep)} ms (λ_max·dt < 1). Errors are relative to the plateau.`);
  log();
  table(['method', 'dt (ms)', 'direction reversals in V_root(t)', 'max error', `error at t = ${ringT} ms`], ringRows);
  out.e3 = { res, hhTraces, sample, ring, lamMax: lam, plateau };
}

// ---------------------------------------------------------------------------
// E4. Forward Euler stability boundary.
// ---------------------------------------------------------------------------
function e4() {
  log('## E4. Explicit (forward Euler) stability boundary');
  log();
  const blowsUp = (model: CableModel, dt: number) => {
    const sim = new Simulator(model, 'fe', dt);
    for (let i = 0; i < model.n; i++) sim.V[i] = i % 2 ? -1 : 1;
    const steps = 4000;
    for (let k = 0; k < steps; k++) {
      sim.step();
      if (!(Math.abs(sim.V[0]!) < 1e3)) return true;
    }
    return Math.max(...sim.V.map(Math.abs)) > 10;
  };
  const rows: (string | number)[][] = [];
  const data: { label: string; dx: number; predicted: number; empirical: number }[] = [];
  const cases: [string, number, SectionSpec[]][] = [];
  for (const nseg of QUICK ? [25, 50, 100] : [25, 50, 100, 200, 400])
    cases.push(['cable 1000 µm', 1000 / nseg, [{ parent: -1, length: 1000, diam: 2, nseg }]]);
  for (const nseg of QUICK ? [4] : [2, 4, 8]) {
    const secs = randomSections(25, nseg, 5);
    cases.push(['random tree (direct)', Math.min(...secs.map((s) => s.length / nseg)), secs]);
  }
  for (const [label, dx, secs] of cases) {
    const model = new CableModel(buildTree(secs, { branch: 'direct' }).morph, PASSIVE);
    const predicted = 2 / model.maxRate(20000);
    let lo = predicted * 0.5, hi = predicted * 1.5;
    for (let k = 0; k < 30; k++) {
      const mid = (lo + hi) / 2;
      if (blowsUp(model, mid)) hi = mid;
      else lo = mid;
    }
    data.push({ label, dx, predicted, empirical: lo });
    rows.push([label, model.n, fmt(dx, 4), fmt(predicted, 5), fmt(lo, 5), fmt(lo / predicted, 5)]);
  }
  log('Largest stable dt from bisection (4000 steps from a checkerboard initial state) vs the prediction 2/λ_max(C⁻¹A):');
  log();
  table(['geometry', 'nodes', 'min dx (µm)', 'predicted dt* (ms)', 'empirical dt* (ms)', 'ratio'], rows);
  const cab = data.filter((d) => d.label.startsWith('cable'));
  log(`On the uniform cable, dt* ∝ dx^${fmt(slope(cab.map((d) => d.dx), cab.map((d) => d.empirical)))}. An implicit step is one O(n) tree solve with no step-size limit.`);
  log();
  out.e4 = data;
}

// ---------------------------------------------------------------------------
// E5. Action-potential propagation: HH 1952's 18.8 m/s and the √d law.
// ---------------------------------------------------------------------------
const AXON: Biophysics = { cm: 1, ra: 35.4, ...HH_LEAK, hh: { ...HH_SQUID, celsius: 18.5 } };
const axonLambdaUm = (diam: number) => lambdaCm(diam, rmFromLeak(HH_LEAK.gLeak), AXON.ra) * 1e4;

/** Conduction velocity (m/s) on a 10λ axon with compartment length ≈ dxUm. */
function velocity(diam: number, dxUm: number, dt: number, method: Method = 'cn'): number {
  const bio = AXON;
  const L = 10 * axonLambdaUm(diam);
  const nseg = Math.round(L / dxUm);
  const sim = new Simulator(new CableModel(cable(L, diam, nseg).morph, bio), method, dt, -65);
  sim.addStim({ node: 0, amp: 2 * Math.pow(diam, 1.5), t0: 0.1, t1: 0.6 }); // ∝ G∞ ∝ d^{3/2}
  const a = Math.round(nseg * 0.3), b = Math.round(nseg * 0.7);
  let ta = NaN, tb = NaN, pa = sim.V[a]!, pb = sim.V[b]!;
  sim.run(12, (s) => {
    const va = s.V[a]!, vb = s.V[b]!;
    if (isNaN(ta) && va >= 0) ta = s.t - (s.dt * va) / (va - pa);
    if (isNaN(tb) && vb >= 0) tb = s.t - (s.dt * vb) / (vb - pb);
    pa = va;
    pb = vb;
  });
  return ((b - a) * (L / nseg) * 1e-6) / ((tb - ta) * 1e-3); // m/s
}

function e5() {
  log('## E5. Action-potential propagation: conduction velocity');
  log();
  const diams = QUICK ? [1, 10, 100, 476] : [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 476, 1000];
  const DX = 10;
  const vs = diams.map((d) => velocity(d, DX, 0.01));
  const vsScaled = diams.map((d) => velocity(d, axonLambdaUm(d) / 40, 0.01));
  log(`Unmyelinated HH axon (squid parameters, 18.5 °C, Ra = 35.4 Ω·cm), 10λ long, CN with dt = 0.01 ms. Velocity is measured between 30% and 70% of the length at the 0 mV crossing.`);
  log(`Grids: a fixed physical dx = ${DX} µm, and one scaled to λ/40.`);
  log();
  table(['diameter (µm)', 'λ (µm)', `v, dx = ${DX} µm (m/s)`, 'v, dx = λ/40 (m/s)', 'v / √d (dx = λ/40)'], diams.map((d, i) => [d, fmt(axonLambdaUm(d), 4), fmt(vs[i]!, 4), fmt(vsScaled[i]!, 4), fmt(vsScaled[i]! / Math.sqrt(d), 5)]));
  const v476 = vsScaled[diams.indexOf(476)]!;
  log(`Fitted exponents: v ∝ d^**${fmt(slope(diams, vs), 4)}** on the fixed grid and d^**${fmt(slope(diams, vsScaled), 6)}** on the λ-scaled grid (theory: 1/2).`);
  log('In the continuous cable equation, v ∝ λ/τ_m with τ_m = R_m C_m independent of d, so √d is exact. The λ-scaled discrete system is exactly self-similar, so it reproduces the exponent to rounding error. That is a consistency check, not evidence. On the fixed grid, small axons are resolved more finely than large ones, and the residual drift measures discretisation error.');
  log(`At the squid giant axon's 476 µm, the model gives **${fmt(v476, 4)} m/s**. Hodgkin & Huxley (1952) computed 18.8 m/s with the same equations (and measured 21.2 m/s).`);
  log();
  const conv: (string | number)[][] = [];
  const convData: { method: string; spl: number; dt: number; v: number }[] = [];
  const grid: [Method, number, number][] = QUICK
    ? [['cn', 10, 0.02], ['cn', 40, 0.01], ['be', 40, 0.01]]
    : [['cn', 10, 0.01], ['cn', 20, 0.01], ['cn', 40, 0.01], ['cn', 80, 0.01], ['cn', 40, 0.04], ['cn', 40, 0.02], ['cn', 40, 0.005], ['be', 40, 0.04], ['be', 40, 0.01], ['be', 40, 0.0025]];
  for (const [method, spl, dt] of grid) {
    const v = velocity(476, axonLambdaUm(476) / spl, dt, method);
    convData.push({ method, spl, dt, v });
    conv.push([method.toUpperCase(), spl, dt, fmt(v, 4)]);
  }
  log('Discretisation sensitivity at d = 476 µm:');
  log();
  table(['method', 'compartments / λ', 'dt (ms)', 'velocity (m/s)'], conv);
  out.e5 = { diams, vs, vsScaled, dx: DX, conv: convData };
}

const t0 = performance.now();
e1();
e2();
e3();
e4();
e5();
log(`_Total runtime ${((performance.now() - t0) / 1000).toFixed(1)} s._`);
mkdirSync('results', { recursive: true });
const suffix = QUICK ? '-quick' : '';
writeFileSync(`results/results${suffix}.json`, JSON.stringify(out));
writeFileSync(`results/results${suffix}.md`, md.join('\n') + '\n');
