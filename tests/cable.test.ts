import { test } from 'node:test';
import assert from 'node:assert/strict';
import { CableModel, PASSIVE, Simulator, type Biophysics, type Method } from '../src/cable.ts';
import { buildTree, cable, rallTree, randomSections, type BranchMode } from '../src/morphology.ts';
import { equivalentCylinder, lambdaCm, rmFromLeak, sealedCableVoltage, treeInputResistance } from '../src/analytic.ts';
import { HH_LEAK, HH_SQUID, steadyState } from '../src/hh.ts';

const RM = rmFromLeak(PASSIVE.gLeak);
const HH_BIO: Biophysics = { cm: 1, ra: 100, ...HH_LEAK, hh: HH_SQUID };

test('sealed-cable steady state matches cosh profile with second-order error', () => {
  const errAt = (nseg: number) => {
    const t = cable(1000, 2, nseg, true);
    const v = new CableModel(t.morph, PASSIVE).steadyState(new Map([[0, 0.05]]));
    let worst = 0;
    for (let k = 0; k < nseg; k++) {
      const exact = sealedCableVoltage((k + 0.5) * (1000 / nseg), 1000, 2, RM, PASSIVE.ra, 0.05);
      worst = Math.max(worst, Math.abs(v[t.first[0]! + k]! - exact) / exact);
    }
    // The zero-area root node sits exactly at x = 0.
    const v0 = sealedCableVoltage(0, 1000, 2, RM, PASSIVE.ra, 0.05);
    return Math.max(worst, Math.abs(v[0]! - v0) / v0);
  };
  const e = [20, 40, 80].map(errAt);
  assert.ok(e[2]! < 1e-3, `error ${e[2]}`);
  for (let i = 1; i < e.length; i++) {
    const ratio = e[i - 1]! / e[i]!;
    assert.ok(ratio > 3.6 && ratio < 4.4, `refinement ratio ${ratio}`);
  }
});

test("a 3/2-rule tree has exactly its equivalent cylinder's input resistance", () => {
  const secs = rallTree(5, 250, 6, 1);
  const a = treeInputResistance(secs, RM, PASSIVE.ra);
  const b = treeInputResistance([equivalentCylinder(secs, RM, PASSIVE.ra)], RM, PASSIVE.ra);
  assert.ok(Math.abs(a - b) / a < 1e-12);
  const broken = secs.map((s, i) => (i ? { ...s, diam: s.diam * 1.2 } : s));
  assert.ok(Math.abs(treeInputResistance(broken, RM, PASSIVE.ra) - a) / a > 0.01);
});

test('junction nodes converge at O(dx²) on branched trees; direct coupling at O(dx)', () => {
  const exact = treeInputResistance(randomSections(15, 1, 9), RM, PASSIVE.ra);
  const err = (mode: BranchMode, nseg: number) => {
    const t = buildTree(randomSections(15, nseg, 9), { branch: mode, rootNode: true });
    return Math.abs(new CableModel(t.morph, PASSIVE).inputResistance(0) - exact) / exact;
  };
  const rj = err('junction', 16) / err('junction', 32);
  const rd = err('direct', 16) / err('direct', 32);
  assert.ok(rj > 3.6 && rj < 4.4, `junction ratio ${rj}`);
  assert.ok(rd > 1.8 && rd < 2.2, `direct ratio ${rd}`);
});

test('backward Euler conserves charge exactly; CN keeps KCL at zero-area nodes', () => {
  const t = buildTree(randomSections(10, 4, 2), { branch: 'junction', rootNode: true });
  const model = new CableModel(t.morph, { ...PASSIVE, eLeak: -70 });
  const dt = 0.05;
  const sim = new Simulator(model, 'be', dt).addStim({ node: 7, amp: 0.3, t0: 0, t1: 5 });
  for (let k = 0; k < 200; k++) {
    const before = Float64Array.from(sim.V);
    const injected = sim.t + dt / 2 < 5 ? 0.3e-3 : 0;
    sim.step();
    let dq = 0, leak = 0;
    for (let i = 0; i < model.n; i++) {
      dq += model.cap[i]! * (sim.V[i]! - before[i]!);
      leak += model.area[i]! * PASSIVE.gLeak * (sim.V[i]! - -70);
    }
    assert.ok(Math.abs(dq - dt * (injected - leak)) < 1e-15, `step ${k}: imbalance ${dq - dt * (injected - leak)}`);
  }

  const cn = new Simulator(model, 'cn', 0.2).addStim({ node: 0, amp: 0.5, t0: 0, t1: 3 });
  const A = model.matrix();
  const y = new Float64Array(model.n);
  cn.run(10, (s) => {
    A.multiply(s.V, y); // axial currents leaving each node
    for (let i = 0; i < model.n; i++) {
      if (model.cap[i] !== 0 || i === 0) continue;
      assert.ok(Math.abs(y[i]!) < 1e-12, `KCL residual ${y[i]} at junction ${i}`);
    }
  });
});

test('temporal order: BE is first order, staggered CN is second order (HH)', () => {
  const t = buildTree(rallTree(2, 150, 3, 8), { branch: 'junction', rootNode: true });
  const probe = t.last[2]!;
  const trace = (method: Method, dt: number) => {
    const sim = new Simulator(new CableModel(t.morph, HH_BIO), method, dt, -65);
    sim.addStim({ node: 0, amp: 0.3, t0: 1, t1: 2 });
    const every = Math.round(0.1 / dt);
    const v: number[] = [];
    let k = 0;
    sim.run(8, (s) => {
      if (++k % every === 0) v.push(s.V[probe]!);
    });
    return v;
  };
  const ref = trace('cn', 0.1 / 128);
  assert.ok(Math.max(...ref) > 0, 'the stimulus should evoke a spike');
  const err = (m: Method, dt: number) => Math.max(...trace(m, dt).map((x, i) => Math.abs(x - ref[i]!)));
  const be = err('be', 0.0125) / err('be', 0.00625);
  const cn = err('cn', 0.025) / err('cn', 0.0125);
  assert.ok(be > 1.7 && be < 2.4, `BE ratio ${be}`);
  assert.ok(cn > 3.4 && cn < 4.6, `CN ratio ${cn}`);
});

test('HH: rest is stationary, a pulse fires a spike, gates stay in [0, 1]', () => {
  const model = new CableModel(cable(100, 20, 3).morph, HH_BIO);
  const rest = new Simulator(model, 'cn', 0.025, -65).run(50);
  for (const v of rest.V) assert.ok(Math.abs(v + 65) < 0.05, `drifted to ${v}`);

  const sim = new Simulator(model, 'cn', 0.025, -65).addStim({ node: 1, amp: 1, t0: 1, t1: 1.5 });
  let peak = -Infinity;
  sim.run(20, (s) => {
    peak = Math.max(peak, s.V[1]!);
    for (const g of [s.m!, s.h!, s.ng!]) for (const x of g) assert.ok(x >= 0 && x <= 1);
  });
  assert.ok(peak > 20 && peak < 50, `peak ${peak}`);
  const ss = steadyState(-65);
  assert.ok(Math.abs(ss.m - 0.0529) < 1e-3 && Math.abs(ss.h - 0.596) < 1e-3 && Math.abs(ss.n - 0.3177) < 1e-3);
});

test('maxRate matches the closed-form top eigenvalue of a uniform sealed chain', () => {
  const n = 60;
  const model = new CableModel(cable(1200, 2, n).morph, PASSIVE);
  const g = model.gAx[1]!, C = model.cap[0]!;
  const exact = (g / C) * (2 + 2 * Math.cos(Math.PI / n)) + PASSIVE.gLeak / PASSIVE.cm;
  assert.ok(Math.abs(model.maxRate(20000) - exact) / exact < 1e-9);
});

test('forward Euler is stable just below 2/λ_max and unstable just above it', () => {
  const model = new CableModel(cable(800, 2, 80).morph, PASSIVE);
  const dtStar = 2 / model.maxRate(20000);
  const grows = (dt: number) => {
    const sim = new Simulator(model, 'fe', dt);
    for (let i = 0; i < model.n; i++) sim.V[i] = i % 2 ? -1 : 1;
    sim.run(dt * 3000);
    return Math.max(...sim.V.map(Math.abs)) > 1;
  };
  assert.equal(grows(0.97 * dtStar), false);
  assert.equal(grows(1.03 * dtStar), true);
  const withJunction = new CableModel(cable(800, 2, 8, true).morph, PASSIVE);
  assert.throws(() => new Simulator(withJunction, 'fe', 1e-3), /forward Euler/);
});

test('squid giant axon conducts near the 18.8 m/s Hodgkin and Huxley computed', () => {
  const bio: Biophysics = { cm: 1, ra: 35.4, ...HH_LEAK, hh: { ...HH_SQUID, celsius: 18.5 } };
  const lam = lambdaCm(476, rmFromLeak(HH_LEAK.gLeak), bio.ra) * 1e4;
  const nseg = 200;
  const L = 10 * lam;
  const sim = new Simulator(new CableModel(cable(L, 476, nseg).morph, bio), 'cn', 0.02, -65);
  sim.addStim({ node: 0, amp: 2 * Math.pow(476, 1.5), t0: 0.1, t1: 0.6 });
  const a = 60, b = 140;
  let ta = NaN, tb = NaN;
  sim.run(12, (s) => {
    if (isNaN(ta) && s.V[a]! >= 0) ta = s.t;
    if (isNaN(tb) && s.V[b]! >= 0) tb = s.t;
  });
  const v = ((b - a) * (L / nseg) * 1e-6) / ((tb - ta) * 1e-3);
  assert.ok(v > 18 && v < 19.2, `velocity ${v} m/s`);
});
