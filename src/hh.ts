// Hodgkin–Huxley (1952) squid-axon kinetics in the modern sign convention
// (resting potential ≈ -65 mV). Rates are in 1/ms at 6.3 °C and are scaled by
// φ = 3^((T - 6.3)/10).

export interface HHParams {
  gNa: number; // mS/cm²
  gK: number; // mS/cm²
  eNa: number; // mV
  eK: number; // mV
  celsius: number;
}

export const HH_SQUID: HHParams = { gNa: 120, gK: 36, eNa: 50, eK: -77, celsius: 6.3 };
/** Leak parameters that go with HH_SQUID. */
export const HH_LEAK = { gLeak: 0.3, eLeak: -54.387 };

/** x / (1 - exp(-x/y)), with its removable singularity at x = 0 filled in. */
function vtrap(x: number, y: number): number {
  const u = x / y;
  return Math.abs(u) < 1e-6 ? y * (1 + u / 2) : x / (1 - Math.exp(-u));
}

export interface Rates {
  am: number; bm: number; ah: number; bh: number; an: number; bn: number;
}

export function rates(v: number, phi: number, out: Rates): Rates {
  out.am = phi * 0.1 * vtrap(v + 40, 10);
  out.bm = phi * 4 * Math.exp(-(v + 65) / 18);
  out.ah = phi * 0.07 * Math.exp(-(v + 65) / 20);
  out.bh = phi / (1 + Math.exp(-(v + 35) / 10));
  out.an = phi * 0.01 * vtrap(v + 55, 10);
  out.bn = phi * 0.125 * Math.exp(-(v + 65) / 80);
  return out;
}

export function phiOf(celsius: number): number {
  return Math.pow(3, (celsius - 6.3) / 10);
}

export function steadyState(v: number): { m: number; h: number; n: number } {
  const r = rates(v, 1, { am: 0, bm: 0, ah: 0, bh: 0, an: 0, bn: 0 });
  return { m: r.am / (r.am + r.bm), h: r.ah / (r.ah + r.bh), n: r.an / (r.an + r.bn) };
}

/**
 * Advance one gate by dt with V held fixed. The gate ODE is linear in that
 * case, so the exponential-Euler update x_inf + (x - x_inf)e^{-dt/τ} is
 * exact. It is unconditionally stable and keeps x in [0, 1].
 */
export function gateStep(x: number, a: number, b: number, dt: number): number {
  const s = a + b;
  const xinf = a / s;
  return xinf + (x - xinf) * Math.exp(-dt * s);
}
