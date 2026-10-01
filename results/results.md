# Results

_Full run._

## E1. Hines solver: correctness, fill-in, and scaling

Max relative difference from dense LU with partial pivoting on 30 backward-Euler cable matrices (random dendrites with junction nodes): **6.43e-16**.

Fill-in (new off-diagonal pairs) from symbolic elimination on uniform random recursive trees:

| n | tree edges | leaves-first (Hines) | root-first | random order | dense n(n-1)/2 |
|---|---|---|---|---|---|
| 100 | 99 | 0 | 2631 | 216 | 4950 |
| 300 | 299 | 0 | 22214 | 912 | 44850 |
| 1000 | 999 | 0 | 250933 | 7923 | 499500 |

Wall-clock time per solve (Node v22.22.2, single thread, median):

| n | Hines (ms) | ns / node |
|---|---|---|
| 1000 | 0.00462 | 4.62 |
| 10000 | 0.0509 | 5.09 |
| 100000 | 0.624 | 6.24 |
| 1000000 | 6.97 | 6.97 |

| n | dense LU (ms) | Hines at same n (ms) | speedup |
|---|---|---|---|
| 100 | 0.243 | 5.42e-4 | 448 |
| 200 | 1.39 | 1.00e-3 | 1.40e+3 |
| 400 | 19.3 | 0.00187 | 1.03e+4 |
| 800 | 141 | 0.00354 | 3.98e+4 |

Fitted exponents: Hines **1.06**, dense LU **3.13**.

## E2. Spatial convergence against exact cable theory

Relative error of the steady-state input resistance at the root, compared with the exact
leaves-up recursion (Rm = 10 kΩ·cm², Ra = 100 Ω·cm). The order is fitted over the four finest grids.

| geometry | branch points | exact R_in (MΩ) | err @ nseg=1 | err @ nseg=128 | observed order |
|---|---|---|---|---|---|
| cable | junction | 253.357 | 0.256 | 1.70e-5 | 2.00 |
| rall | junction | 95.4562 | 0.0128 | 7.87e-7 | 2.00 |
| rall | direct | 95.4562 | 0.0311 | 3.32e-4 | 0.993 |
| random | junction | 83.5946 | 0.0112 | 6.92e-7 | 2.00 |
| random | direct | 83.5946 | 0.0485 | 4.46e-4 | 0.996 |

**Rall's equivalent cylinder.** A 4-generation binary tree obeying d_p^{3/2} = Σ d_c^{3/2}, with equal electrotonic path lengths, has exact R_in = 95.456165 MΩ. Its equivalent cylinder (L = 1200.0 µm) has 95.456165 MΩ.
The transient check: under a 0.1 nA step, compare the simulated root voltage of the tree and the cylinder (BE, dt = 0.002 ms, 30 ms), as max |ΔV| / plateau:

| compartments in root section | tree vs cylinder | 3/2 rule broken (daughters ×1.3 diameter) |
|---|---|---|
| 4 | 6.89% | 16.4% |
| 8 | 2.09% | 15.6% |
| 16 | 0.488% | 15.4% |
| 32 | 0.117% | 15.4% |

With the 3/2 rule, the tree–cylinder gap is a discretisation artifact that shrinks under refinement. When the rule is broken, a real, grid-independent difference remains.

## E3. Temporal convergence: backward Euler vs staggered Crank–Nicolson

Max |V − V_ref| (mV) at a distal tip of a 3-generation tree, sampled every 0.1 ms over 15 ms (2 ms current step).
Reference: CN at dt = 1.95e-4 ms. HH is the squid axon at 6.3 °C on every compartment.

| membrane | method | err @ dt=0.1 | err @ dt=0.00625 | observed order |
|---|---|---|---|---|
| passive | BE | 0.0916 | 0.00638 | 0.962 |
| passive | CN | 0.00387 | 1.66e-5 | 1.98 |
| HH | BE | 61.2 | 4.54 | 0.952 |
| HH | CN | 13.2 | 0.0525 | 1.99 |

**Crank–Nicolson ringing.** Root voltage of a sealed 1000 µm cable (nseg = 200) under a 0.1 nA step. The exact response rises monotonically.
The stiffest mode has λ_max ≈ 8.00e+3 /ms. At dt = 0.5 ms its CN amplification factor is (1 − λdt/2)/(1 + λdt/2) ≈ -0.9990, so it flips sign every step and barely decays.
The reference is CN at dt = 1.25e-4 ms (λ_max·dt < 1). Errors are relative to the plateau.

| method | dt (ms) | direction reversals in V_root(t) | max error | error at t = 10 ms |
|---|---|---|---|---|
| BE | 0.500 | 0 | 2.35% | 0.313% |
| CN | 0.500 | 19 | 6.02% | 0.956% |
| BE | 0.0500 | 0 | 0.493% | 0.251% |
| CN | 0.0500 | 199 | 2.13% | 0.00259% |

## E4. Explicit (forward Euler) stability boundary

Largest stable dt from bisection (4000 steps from a checkerboard initial state) vs the prediction 2/λ_max(C⁻¹A):

| geometry | nodes | min dx (µm) | predicted dt* (ms) | empirical dt* (ms) | ratio |
|---|---|---|---|---|---|
| cable 1000 µm | 25 | 40 | 0.016050 | 0.016055 | 1.0003 |
| cable 1000 µm | 50 | 20 | 0.0040031 | 0.0040042 | 1.0003 |
| cable 1000 µm | 100 | 10 | 0.0010002 | 0.0010005 | 1.0003 |
| cable 1000 µm | 200 | 5 | 2.5001e-4 | 2.5008e-4 | 1.0003 |
| cable 1000 µm | 400 | 2.500 | 6.2501e-5 | 6.2518e-5 | 1.0003 |
| random tree (direct) | 50 | 27.38 | 0.012859 | 0.012863 | 1.0003 |
| random tree (direct) | 100 | 13.69 | 0.0021764 | 0.0021770 | 1.0003 |
| random tree (direct) | 200 | 6.844 | 4.8986e-4 | 4.8998e-4 | 1.0003 |

On the uniform cable, dt* ∝ dx^2.00. An implicit step is one O(n) tree solve with no step-size limit.

## E5. Action-potential propagation: conduction velocity

Unmyelinated HH axon (squid parameters, 18.5 °C, Ra = 35.4 Ω·cm), 10λ long, CN with dt = 0.01 ms. Velocity is measured between 30% and 70% of the length at the 0 mV crossing.
Grids: a fixed physical dx = 10 µm, and one scaled to λ/40.

| diameter (µm) | λ (µm) | v, dx = 10 µm (m/s) | v, dx = λ/40 (m/s) | v / √d (dx = λ/40) |
|---|---|---|---|---|
| 0.500 | 343.1 | 0.6059 | 0.6060 | 0.85707 |
| 1 | 485.2 | 0.8572 | 0.8571 | 0.85707 |
| 2 | 686.2 | 1.213 | 1.212 | 0.85707 |
| 5 | 1085 | 1.917 | 1.916 | 0.85707 |
| 10 | 1534 | 2.712 | 2.710 | 0.85707 |
| 20 | 2170 | 3.835 | 3.833 | 0.85707 |
| 50 | 3431 | 6.064 | 6.060 | 0.85707 |
| 100 | 4852 | 8.576 | 8.571 | 0.85707 |
| 200 | 6862 | 12.13 | 12.12 | 0.85707 |
| 476 | 1.059e+4 | 18.71 | 18.70 | 0.85707 |
| 1000 | 1.534e+4 | 27.12 | 27.10 | 0.85707 |

Fitted exponents: v ∝ d^**0.5001** on the fixed grid and d^**0.500000** on the λ-scaled grid (theory: 1/2).
In the continuous cable equation, v ∝ λ/τ_m with τ_m = R_m C_m independent of d, so √d is exact. The λ-scaled discrete system is exactly self-similar, so it reproduces the exponent to rounding error. That is a consistency check, not evidence. On the fixed grid, small axons are resolved more finely than large ones, and the residual drift measures discretisation error.
At the squid giant axon's 476 µm, the model gives **18.70 m/s**. Hodgkin & Huxley (1952) computed 18.8 m/s with the same equations (and measured 21.2 m/s).

Discretisation sensitivity at d = 476 µm:

| method | compartments / λ | dt (ms) | velocity (m/s) |
|---|---|---|---|
| CN | 10 | 0.0100 | 18.54 |
| CN | 20 | 0.0100 | 18.67 |
| CN | 40 | 0.0100 | 18.70 |
| CN | 80 | 0.0100 | 18.71 |
| CN | 40 | 0.0400 | 18.36 |
| CN | 40 | 0.0200 | 18.63 |
| CN | 40 | 0.00500 | 18.72 |
| BE | 40 | 0.0400 | 18.00 |
| BE | 40 | 0.0100 | 18.58 |
| BE | 40 | 0.00250 | 18.69 |

_Total runtime 11.3 s._
