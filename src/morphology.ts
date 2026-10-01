// Neuron morphologies as rooted trees of nodes in Hines order (parent[i] < i).
//
// A node is either a compartment (a cylinder of length L and diameter d, with
// membrane on its side surface) or a zero-area junction node (length 0). Its
// voltage is defined at the compartment centre or at the junction point. The
// axial path from node i to its parent runs through half of each cylinder, so
// a zero-length node adds no resistance of its own.

export interface Morphology {
  readonly n: number;
  readonly parent: Int32Array; // parent[0] = -1 and parent[i] < i for i > 0
  readonly length: Float64Array; // µm; 0 marks a zero-area junction node
  readonly diam: Float64Array; // µm
}

export interface SectionSpec {
  parent: number; // index of the parent section (-1 for the root); must be < own index
  length: number; // µm
  diam: number; // µm
  nseg: number; // number of compartments in the section
}

/**
 * How children attach at a branch point:
 *  - 'direct':   each child's first compartment couples to the parent's last
 *                compartment through (R_parent/2 + R_child/2). With k children,
 *                the parent's distal half-resistance appears k times in parallel
 *                when it should appear once in series, which gives an O(dx) error.
 *  - 'junction': a zero-area node sits at the branch point, so each half-
 *                resistance appears exactly once (NEURON's approach). O(dx²).
 */
export type BranchMode = 'direct' | 'junction';

export interface BuildOptions {
  branch?: BranchMode;
  /** Put a zero-area node at the proximal end of the root section, where x = 0. */
  rootNode?: boolean;
}

export interface BuiltTree {
  morph: Morphology;
  /** First and last compartment node of each section. */
  first: Int32Array;
  last: Int32Array;
  /** Node index at section-relative position x in [0,1] (nearest compartment centre). */
  nodeAt(section: number, x: number): number;
}

export function buildTree(sections: readonly SectionSpec[], opts: BuildOptions = {}): BuiltTree {
  const branch = opts.branch ?? 'junction';
  const hasChildren = new Array<boolean>(sections.length).fill(false);
  sections.forEach((s, i) => {
    if (i === 0 ? s.parent !== -1 : s.parent < 0 || s.parent >= i)
      throw new Error(`section ${i}: parent must precede it (got ${s.parent})`);
    if (!(s.nseg >= 1 && s.length > 0 && s.diam > 0)) throw new Error(`section ${i}: bad geometry`);
    if (s.parent >= 0) hasChildren[s.parent] = true;
  });

  const parent: number[] = [];
  const length: number[] = [];
  const diam: number[] = [];
  const add = (p: number, L: number, d: number): number => {
    parent.push(p);
    length.push(L);
    diam.push(d);
    return parent.length - 1;
  };

  const attach = new Int32Array(sections.length);
  const first = new Int32Array(sections.length);
  const last = new Int32Array(sections.length);
  sections.forEach((s, i) => {
    let p: number;
    if (s.parent === -1) p = opts.rootNode ? add(-1, 0, s.diam) : -1;
    else p = attach[s.parent]!;
    const dx = s.length / s.nseg;
    for (let k = 0; k < s.nseg; k++) {
      p = add(p, dx, s.diam);
      if (k === 0) first[i] = p;
    }
    last[i] = p;
    attach[i] = branch === 'junction' && hasChildren[i] ? add(p, 0, s.diam) : p;
  });

  const morph: Morphology = {
    n: parent.length,
    parent: Int32Array.from(parent),
    length: Float64Array.from(length),
    diam: Float64Array.from(diam),
  };
  return {
    morph,
    first,
    last,
    nodeAt(section, x) {
      const s = sections[section]!;
      const k = Math.min(s.nseg - 1, Math.max(0, Math.floor(x * s.nseg)));
      return first[section]! + k;
    },
  };
}

/** An unbranched cylinder. */
export function cable(lengthUm: number, diamUm: number, nseg: number, rootNode = false): BuiltTree {
  return buildTree([{ parent: -1, length: lengthUm, diam: diamUm, nseg }], { rootNode });
}

/** Mulberry32: small deterministic PRNG. */
export function rng(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Random node-level tree (uniform random recursive tree) for solver tests. */
export function randomNodeTree(n: number, seed: number): Morphology {
  const r = rng(seed);
  const parent = new Int32Array(n);
  const length = new Float64Array(n);
  const diam = new Float64Array(n);
  parent[0] = -1;
  for (let i = 0; i < n; i++) {
    if (i > 0) parent[i] = Math.floor(r() * i);
    length[i] = 5 + 45 * r();
    diam[i] = 0.5 + 4.5 * r();
  }
  return { n, parent, length, diam };
}

/**
 * Random dendritic tree of nSections sections. Each section after the first
 * hangs off the distal end of an earlier section. Diameters taper with depth,
 * and the result deliberately does not satisfy Rall's 3/2 rule.
 */
export function randomSections(nSections: number, nseg: number, seed: number): SectionSpec[] {
  const r = rng(seed);
  const depth: number[] = [0];
  const out: SectionSpec[] = [{ parent: -1, length: 100 + 200 * r(), diam: 3 + 2 * r(), nseg }];
  for (let i = 1; i < nSections; i++) {
    const p = Math.floor(r() * i);
    depth.push(depth[p]! + 1);
    out.push({
      parent: p,
      length: 50 + 250 * r(),
      diam: Math.max(0.3, out[p]!.diam * (0.5 + 0.4 * r())),
      nseg,
    });
  }
  return out;
}

/**
 * Symmetric binary tree that satisfies Rall's 3/2 rule: at every branch point
 * d_parent^{3/2} = 2 d_child^{3/2}. Child lengths are chosen so each generation
 * has the same electrotonic length (L ∝ sqrt(d)), so all terminals end at the
 * same electrotonic distance and the tree collapses to an equivalent cylinder.
 */
export function rallTree(levels: number, rootLength: number, rootDiam: number, nseg: number): SectionSpec[] {
  const out: SectionSpec[] = [{ parent: -1, length: rootLength, diam: rootDiam, nseg }];
  let frontier = [0];
  for (let lvl = 1; lvl < levels; lvl++) {
    const next: number[] = [];
    for (const p of frontier) {
      const d = out[p]!.diam * Math.pow(2, -2 / 3);
      const L = out[p]!.length * Math.sqrt(d / out[p]!.diam);
      for (let c = 0; c < 2; c++) {
        out.push({ parent: p, length: L, diam: d, nseg });
        next.push(out.length - 1);
      }
    }
    frontier = next;
  }
  return out;
}
