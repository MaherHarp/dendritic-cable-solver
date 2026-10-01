// Closed-form passive cable theory (Rall), used as an independent reference
// for the discretised model.
//
// For a cylinder of diameter d with specific membrane resistance Rm (Ω·cm²)
// and axial resistivity Ra (Ω·cm):
//   λ   = sqrt((d/4) Rm / Ra)              length constant
//   G∞  = 1 / (r_a λ),  r_a = 4 Ra/(π d²)  input conductance of a semi-infinite cable
// A finite cylinder of electrotonic length X = L/λ that ends in a load
// conductance G_L has input conductance
//   G_in = G∞ (G_L/G∞ + tanh X) / (1 + (G_L/G∞) tanh X).
// Applying this from the leaves up (children in parallel form the parent's
// load) gives the exact steady-state input resistance of any passive tree.

import type { SectionSpec } from './morphology.ts';

/** Specific membrane resistance (Ω·cm²) for a leak conductance in mS/cm². */
export function rmFromLeak(gLeakMs: number): number {
  return 1000 / gLeakMs;
}

/** Length constant in cm. */
export function lambdaCm(diamUm: number, rm: number, ra: number): number {
  return Math.sqrt(((diamUm * 1e-4) / 4) * (rm / ra));
}

/** Semi-infinite input conductance in S. */
export function gInf(diamUm: number, rm: number, ra: number): number {
  const d = diamUm * 1e-4;
  const raPerCm = (4 * ra) / (Math.PI * d * d);
  return 1 / (raPerCm * lambdaCm(diamUm, rm, ra));
}

export function cylinderInputConductance(gLoad: number, gInfS: number, X: number): number {
  const t = Math.tanh(X);
  const B = gLoad / gInfS;
  return (gInfS * (B + t)) / (1 + B * t);
}

/** Exact input conductance (S) seen at the proximal end of each section. */
export function sectionInputConductances(sections: readonly SectionSpec[], rm: number, ra: number): Float64Array {
  const gin = new Float64Array(sections.length);
  const load = new Float64Array(sections.length);
  for (let i = sections.length - 1; i >= 0; i--) {
    const s = sections[i]!;
    const X = (s.length * 1e-4) / lambdaCm(s.diam, rm, ra);
    gin[i] = cylinderInputConductance(load[i]!, gInf(s.diam, rm, ra), X);
    if (s.parent >= 0) load[s.parent] += gin[i]!;
  }
  return gin;
}

/** Exact input resistance at the proximal end of the root section, in MΩ. */
export function treeInputResistance(sections: readonly SectionSpec[], rm: number, ra: number): number {
  return 1 / sectionInputConductances(sections, rm, ra)[0]! / 1e6;
}

/**
 * Steady-state voltage (relative to rest, mV) at distance x (µm) along a sealed
 * finite cylinder of length L (µm) when I (nA) is injected at x = 0:
 *   V(X) = I R∞ cosh(L_e - X) / sinh(L_e),  with R∞ = 1/G∞.
 */
export function sealedCableVoltage(xUm: number, lengthUm: number, diamUm: number, rm: number, ra: number, iNa: number): number {
  const lam = lambdaCm(diamUm, rm, ra);
  const Le = (lengthUm * 1e-4) / lam;
  const X = (xUm * 1e-4) / lam;
  const rInfMOhm = 1 / gInf(diamUm, rm, ra) / 1e6;
  return (iNa * rInfMOhm * Math.cosh(Le - X)) / Math.sinh(Le);
}

/**
 * Rall's equivalent cylinder for a tree that satisfies the 3/2 rule with equal
 * electrotonic path lengths. It has the root's diameter and the tree's total
 * electrotonic length, expressed as a physical length at the root's λ.
 */
export function equivalentCylinder(sections: readonly SectionSpec[], rm: number, ra: number): SectionSpec {
  const root = sections[0]!;
  let X = 0;
  let s: number | undefined = 0;
  const childOf = (p: number) => sections.findIndex((c) => c.parent === p);
  while (s !== undefined && s >= 0) {
    const sec = sections[s]!;
    X += (sec.length * 1e-4) / lambdaCm(sec.diam, rm, ra);
    const c = childOf(s);
    s = c >= 0 ? c : undefined;
  }
  return { parent: -1, length: (X * lambdaCm(root.diam, rm, ra)) / 1e-4, diam: root.diam, nseg: root.nseg };
}
