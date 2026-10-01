// The Hines algorithm (Hines 1984): Gaussian elimination for a matrix whose
// sparsity pattern is a tree, in O(n) time and with no fill-in.
//
// Row i has a diagonal d[i] and at most one "parent" coupling per direction:
//   lower[i] = A[i][parent[i]]   (row i, parent column)
//   upper[i] = A[parent[i]][i]   (parent row, column i)
// plus the couplings to its children, which are stored as those children's
// upper entries. If parent[i] < i for every i, then the nodes n-1, ..., 1 are
// always leaves of the not-yet-eliminated subtree when they are eliminated.
// Eliminating a leaf changes only its parent's diagonal and right-hand side,
// so no new nonzeros appear. Back-substitution then runs root-to-leaves.

export class TreeMatrix {
  readonly n: number;
  readonly parent: Int32Array;
  readonly diag: Float64Array;
  readonly lower: Float64Array;
  readonly upper: Float64Array;
  private readonly d: Float64Array; // scratch copy of diag, destroyed by solve()

  constructor(parent: Int32Array) {
    const n = parent.length;
    if (n === 0 || parent[0] !== -1) throw new Error('node 0 must be the root');
    for (let i = 1; i < n; i++) {
      const p = parent[i]!;
      if (p < 0 || p >= i) throw new Error(`not in Hines order: parent[${i}] = ${p}`);
    }
    this.n = n;
    this.parent = parent;
    this.diag = new Float64Array(n);
    this.lower = new Float64Array(n);
    this.upper = new Float64Array(n);
    this.d = new Float64Array(n);
  }

  /** y = A x (used for explicit stepping and for residual checks). */
  multiply(x: Float64Array, y: Float64Array): Float64Array {
    const { n, parent, diag, lower, upper } = this;
    for (let i = 0; i < n; i++) y[i] = diag[i]! * x[i]!;
    for (let i = 1; i < n; i++) {
      const p = parent[i]!;
      y[i] += lower[i]! * x[p]!;
      y[p] += upper[i]! * x[i]!;
    }
    return y;
  }

  /**
   * Solve A x = b in place: b is overwritten with x. This takes 2(n-1)
   * multiply-adds plus n divisions in each sweep, and no pivoting. That is safe
   * for the diagonally dominant M-matrices produced by cable discretisations.
   */
  solve(b: Float64Array): Float64Array {
    const { n, parent, lower, upper } = this;
    const d = this.d;
    d.set(this.diag);
    for (let i = n - 1; i > 0; i--) {
      const p = parent[i]!;
      const f = upper[i]! / d[i]!;
      d[p] -= f * lower[i]!;
      b[p] -= f * b[i]!;
    }
    b[0] /= d[0]!;
    for (let i = 1; i < n; i++) b[i] = (b[i]! - lower[i]! * b[parent[i]!]!) / d[i]!;
    return b;
  }

  /** Dense copy, for verification against general-purpose LU. */
  toDense(): Float64Array[] {
    const rows = Array.from({ length: this.n }, () => new Float64Array(this.n));
    for (let i = 0; i < this.n; i++) rows[i]![i] = this.diag[i]!;
    for (let i = 1; i < this.n; i++) {
      const p = this.parent[i]!;
      rows[i]![p] = this.lower[i]!;
      rows[p]![i] = this.upper[i]!;
    }
    return rows;
  }
}
