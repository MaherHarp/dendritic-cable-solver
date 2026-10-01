// Reference solvers used to check the Hines solver and to show what it avoids:
// dense LU with partial pivoting (O(n³)) and symbolic elimination, which counts
// the fill-in a given elimination order would create on a sparse graph.

/** Solve A x = b by LU with partial pivoting. A and b are not modified. */
export function denseSolve(A: readonly Float64Array[], b: Float64Array): Float64Array {
  const n = A.length;
  const M = A.map((r) => Float64Array.from(r));
  const x = Float64Array.from(b);
  for (let k = 0; k < n; k++) {
    let piv = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i]![k]!) > Math.abs(M[piv]![k]!)) piv = i;
    if (M[piv]![k] === 0) throw new Error('singular matrix');
    if (piv !== k) {
      [M[k], M[piv]] = [M[piv]!, M[k]!];
      [x[k], x[piv]] = [x[piv]!, x[k]!];
    }
    const rk = M[k]!;
    const akk = rk[k]!;
    for (let i = k + 1; i < n; i++) {
      const ri = M[i]!;
      const f = ri[k]! / akk;
      if (f === 0) continue;
      for (let j = k; j < n; j++) ri[j] -= f * rk[j]!;
      x[i] -= f * x[k]!;
    }
  }
  for (let i = n - 1; i >= 0; i--) {
    const ri = M[i]!;
    let s = x[i]!;
    for (let j = i + 1; j < n; j++) s -= ri[j]! * x[j]!;
    x[i] = s / ri[i]!;
  }
  return x;
}

/**
 * Symbolic Gaussian elimination on an undirected graph. Eliminating vertex v
 * makes its remaining neighbours a clique, and every edge added this way is
 * one fill-in nonzero (counted once per symmetric pair). `order` lists the
 * vertices in elimination order.
 */
export function countFill(parent: Int32Array, order: ArrayLike<number>): number {
  const n = parent.length;
  const adj = Array.from({ length: n }, () => new Set<number>());
  for (let i = 0; i < n; i++) {
    const p = parent[i]!;
    if (p >= 0) {
      adj[i]!.add(p);
      adj[p]!.add(i);
    }
  }
  const gone = new Uint8Array(n);
  let fill = 0;
  for (let k = 0; k < order.length; k++) {
    const v = order[k]!;
    const nb = [...adj[v]!].filter((u) => !gone[u]);
    for (let a = 0; a < nb.length; a++) {
      const sa = adj[nb[a]!]!;
      for (let b = a + 1; b < nb.length; b++) {
        if (!sa.has(nb[b]!)) {
          sa.add(nb[b]!);
          adj[nb[b]!]!.add(nb[a]!);
          fill++;
        }
      }
    }
    gone[v] = 1;
  }
  return fill;
}

/** Elimination orders worth comparing on a tree in Hines order. */
export function orders(n: number, rand: () => number): Record<'hines' | 'rootFirst' | 'random', Int32Array> {
  const hines = new Int32Array(n);
  const rootFirst = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    hines[i] = n - 1 - i; // leaves first: the order TreeMatrix.solve uses
    rootFirst[i] = i;
  }
  const random = Int32Array.from(rootFirst);
  for (let i = n - 1; i > 0; i--) {
    const j = Math.floor(rand() * (i + 1));
    [random[i], random[j]] = [random[j]!, random[i]!];
  }
  return { hines, rootFirst, random };
}
