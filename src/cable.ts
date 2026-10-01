// Compartmental cable model and time integrators.
//
// Units: V mV, t ms, conductance mS, capacitance µF, current µA (stimuli are
// given in nA). These are consistent: mS·mV = µA and µF·mV/ms = µA. Each node
// obeys Kirchhoff's current law in absolute (not per-area) form:
//
//   C_i dV_i/dt = Σ_j g_ij (V_j - V_i) - Σ_k g_k,i (V_i - E_k) + I_i(t)
//
// This form makes the axial part a symmetric tree-structured M-matrix. Given
// the gate values, every membrane current is linear in V, so an implicit
// step is a single tree solve with no Newton iterations.

import { TreeMatrix } from './hines.ts';
import type { Morphology } from './morphology.ts';
import { type HHParams, type Rates, gateStep, phiOf, rates, steadyState as hhSteady } from './hh.ts';

export interface Biophysics {
  cm: number; // µF/cm²
  ra: number; // Ω·cm
  gLeak: number; // mS/cm²
  eLeak: number; // mV
  hh?: HHParams | undefined; // Na/K channels on every compartment with membrane
}

export const PASSIVE: Biophysics = { cm: 1, ra: 100, gLeak: 0.1, eLeak: 0 };

export class CableModel {
  readonly n: number;
  readonly area: Float64Array; // cm²
  readonly cap: Float64Array; // µF
  readonly gAx: Float64Array; // mS, from node i to parent[i]; gAx[0] = 0
  readonly diagAx: Float64Array; // Σ of axial conductances incident on node i
  readonly morph: Morphology;
  readonly bio: Biophysics;

  constructor(morph: Morphology, bio: Biophysics) {
    this.morph = morph;
    this.bio = bio;
    const { n, parent, length, diam } = morph;
    this.n = n;
    this.area = new Float64Array(n);
    this.cap = new Float64Array(n);
    this.gAx = new Float64Array(n);
    this.diagAx = new Float64Array(n);
    // Half-cylinder axial resistance in kΩ: Ra (L/2) / (π r²), with µm → cm.
    const halfR = (i: number) => {
      const r = diam[i]! * 0.5e-4;
      return (bio.ra * length[i]! * 0.5e-4) / (Math.PI * r * r) / 1000;
    };
    for (let i = 0; i < n; i++) {
      this.area[i] = Math.PI * diam[i]! * length[i]! * 1e-8;
      this.cap[i] = bio.cm * this.area[i]!;
    }
    for (let i = 1; i < n; i++) {
      const R = halfR(i) + halfR(parent[i]!);
      if (!(R > 0)) throw new Error(`zero axial resistance between ${i} and its parent`);
      this.gAx[i] = 1 / R;
      this.diagAx[i] += this.gAx[i]!;
      this.diagAx[parent[i]!] += this.gAx[i]!;
    }
  }

  /** Tree matrix with the axial couplings filled in; diag set to diagAx + extra. */
  matrix(extraDiag?: Float64Array): TreeMatrix {
    const A = new TreeMatrix(this.morph.parent);
    for (let i = 1; i < this.n; i++) A.lower[i] = A.upper[i] = -this.gAx[i]!;
    for (let i = 0; i < this.n; i++) A.diag[i] = this.diagAx[i]! + (extraDiag ? extraDiag[i]! : 0);
    return A;
  }

  /** Passive steady state (dt → ∞) under constant injected currents in nA. */
  steadyState(injectNa: ReadonlyMap<number, number>): Float64Array {
    const gl = this.area.map((a) => a * this.bio.gLeak);
    const A = this.matrix(gl);
    const b = gl.map((g) => g * this.bio.eLeak);
    for (const [node, amp] of injectNa) b[node] += amp * 1e-3;
    return A.solve(b);
  }

  /** Input resistance at `node` in MΩ (mV per nA) from the passive steady state. */
  inputResistance(node: number): number {
    const v = this.steadyState(new Map([[node, 1]]));
    return v[node]! - this.bio.eLeak;
  }

  /**
   * Largest eigenvalue of C⁻¹(L + G_leak) in 1/ms, by power iteration on the
   * similar symmetric matrix C^{-1/2}(L + G)C^{-1/2}. Forward Euler is stable
   * exactly when dt < 2 / λ_max.
   */
  maxRate(iters = 3000): number {
    if (this.cap.some((c) => c === 0)) throw new Error('zero-capacitance nodes have infinite rate');
    const A = this.matrix(this.area.map((a) => a * this.bio.gLeak));
    const s = this.cap.map((c) => 1 / Math.sqrt(c));
    const n = this.n;
    // Start near the most oscillatory mode, which is the top mode of a chain.
    const x = Float64Array.from({ length: n }, (_, i) => (i % 2 ? -1 : 1) * (1 + 0.1 * Math.sin(i)));
    const z = new Float64Array(n);
    const y = new Float64Array(n);
    let lam = 0;
    for (let k = 0; k < iters; k++) {
      let xx = 0;
      for (let i = 0; i < n; i++) xx += x[i]! * x[i]!;
      const inv = 1 / Math.sqrt(xx);
      for (let i = 0; i < n; i++) z[i] = s[i]! * x[i]! * inv;
      A.multiply(z, y);
      lam = 0; // Rayleigh quotient of the normalised x
      for (let i = 0; i < n; i++) {
        y[i] *= s[i]!;
        lam += y[i]! * x[i]! * inv;
      }
      x.set(y);
    }
    return lam;
  }
}

/**
 *  - 'be': backward Euler, first order, L-stable.
 *  - 'cn': Crank–Nicolson as NEURON implements it: a backward-Euler half step
 *          to t + dt/2, then V(t+dt) = 2 V(t+dt/2) - V(t). Gates are staggered
 *          at half steps, which makes the scheme second order. It is
 *          A-stable but not L-stable, so stiff modes ring at large dt.
 *  - 'fe': forward Euler, conditionally stable (dt < 2/λ_max).
 */
export type Method = 'be' | 'cn' | 'fe';

export interface Stim {
  node: number;
  amp: number; // nA
  t0: number; // ms, inclusive
  t1: number; // ms, exclusive
}

export class Simulator {
  t = 0;
  readonly V: Float64Array;
  readonly m: Float64Array | null = null;
  readonly h: Float64Array | null = null;
  readonly ng: Float64Array | null = null;
  readonly stims: Stim[] = [];
  private readonly A: TreeMatrix;
  private readonly b: Float64Array;
  private readonly G: Float64Array; // total membrane conductance per node, mS
  private readonly GE: Float64Array; // Σ g_k E_k per node, µA
  private readonly I: Float64Array;
  private readonly phi: number;
  private readonly r: Rates = { am: 0, bm: 0, ah: 0, bh: 0, an: 0, bn: 0 };

  readonly model: CableModel;
  readonly method: Method;
  readonly dt: number;

  constructor(model: CableModel, method: Method, dt: number, v0: number = model.bio.eLeak) {
    this.model = model;
    this.method = method;
    this.dt = dt;
    const n = model.n;
    if (method === 'fe' && model.cap.some((c) => c === 0))
      throw new Error('forward Euler needs capacitance on every node (no junction nodes)');
    this.V = new Float64Array(n).fill(v0);
    this.A = model.matrix();
    this.b = new Float64Array(n);
    this.G = new Float64Array(n);
    this.GE = new Float64Array(n);
    this.I = new Float64Array(n);
    this.phi = phiOf(model.bio.hh?.celsius ?? 6.3);
    if (model.bio.hh) {
      const ss = hhSteady(v0);
      this.m = new Float64Array(n).fill(ss.m);
      this.h = new Float64Array(n).fill(ss.h);
      this.ng = new Float64Array(n).fill(ss.n);
    }
  }

  addStim(s: Stim): this {
    this.stims.push(s);
    return this;
  }

  /** Total membrane + stimulus current balance is evaluated at the step midpoint. */
  private currents(): void {
    const { model, G, GE, I, V } = this;
    const { area } = model;
    const { gLeak, eLeak, hh } = model.bio;
    I.fill(0);
    const tm = this.t + this.dt / 2;
    for (const s of this.stims) if (tm >= s.t0 && tm < s.t1) I[s.node] += s.amp * 1e-3;
    if (hh) {
      const m = this.m!, h = this.h!, ng = this.ng!, r = this.r, dt = this.dt, phi = this.phi;
      for (let i = 0; i < model.n; i++) {
        rates(V[i]!, phi, r);
        m[i] = gateStep(m[i]!, r.am, r.bm, dt);
        h[i] = gateStep(h[i]!, r.ah, r.bh, dt);
        ng[i] = gateStep(ng[i]!, r.an, r.bn, dt);
      }
      for (let i = 0; i < model.n; i++) {
        const a = area[i]!;
        const gna = hh.gNa * m[i]! ** 3 * h[i]!;
        const gk = hh.gK * ng[i]! ** 4;
        G[i] = a * (gLeak + gna + gk);
        GE[i] = a * (gLeak * eLeak + gna * hh.eNa + gk * hh.eK);
      }
    } else {
      for (let i = 0; i < model.n; i++) {
        G[i] = area[i]! * gLeak;
        GE[i] = G[i]! * eLeak;
      }
    }
  }

  step(): void {
    this.currents();
    const { model, A, b, G, GE, I, V, dt } = this;
    const { cap, diagAx, n } = model;
    if (this.method === 'fe') {
      for (let i = 0; i < n; i++) A.diag[i] = diagAx[i]! + G[i]!;
      A.multiply(V, b);
      for (let i = 0; i < n; i++) V[i] += (dt / cap[i]!) * (GE[i]! + I[i]! - b[i]!);
    } else {
      const c = this.method === 'be' ? 1 : 2;
      for (let i = 0; i < n; i++) {
        const cdt = (c * cap[i]!) / dt;
        A.diag[i] = diagAx[i]! + G[i]! + cdt;
        b[i] = cdt * V[i]! + GE[i]! + I[i]!;
      }
      A.solve(b);
      if (c === 1) V.set(b);
      else for (let i = 0; i < n; i++) V[i] = 2 * b[i]! - V[i]!;
    }
    this.t += dt;
  }

  /** Step until t reaches tEnd (rounded to whole steps), calling `observe` after each step. */
  run(tEnd: number, observe?: (sim: this) => void): this {
    const steps = Math.round((tEnd - this.t) / this.dt);
    for (let k = 0; k < steps; k++) {
      this.step();
      observe?.(this);
    }
    return this;
  }
}
