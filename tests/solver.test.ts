import { test } from 'node:test';
import assert from 'node:assert/strict';
import { TreeMatrix } from '../src/hines.ts';
import { denseSolve, countFill, orders } from '../src/dense.ts';
import { randomNodeTree, rng, buildTree, randomSections } from '../src/morphology.ts';

function randomSystem(n: number, seed: number, symmetric: boolean): { A: TreeMatrix; b: Float64Array } {
  const r = rng(seed);
  const m = randomNodeTree(n, seed);
  const A = new TreeMatrix(m.parent);
  for (let i = 1; i < n; i++) {
    A.lower[i] = -(0.1 + r());
    A.upper[i] = symmetric ? A.lower[i]! : -(0.1 + r());
  }
  // Strict diagonal dominance in columns and rows, like a cable matrix with C/dt > 0.
  for (let i = 0; i < n; i++) A.diag[i] = 0.05 + r();
  for (let i = 1; i < n; i++) {
    const p = m.parent[i]!;
    const s = Math.abs(A.lower[i]!) + Math.abs(A.upper[i]!);
    A.diag[i] += s;
    A.diag[p] += s;
  }
  const b = Float64Array.from({ length: n }, () => r() * 2 - 1);
  return { A, b };
}

test('Hines solve matches dense LU on random trees (symmetric and not)', () => {
  for (const symmetric of [true, false]) {
    for (let seed = 1; seed <= 20; seed++) {
      const n = 2 + ((seed * 37) % 150);
      const { A, b } = randomSystem(n, seed, symmetric);
      const ref = denseSolve(A.toDense(), b);
      const x = A.solve(Float64Array.from(b));
      let err = 0;
      let scale = 0;
      for (let i = 0; i < n; i++) {
        err = Math.max(err, Math.abs(x[i]! - ref[i]!));
        scale = Math.max(scale, Math.abs(ref[i]!));
      }
      assert.ok(err / scale < 1e-12, `seed ${seed} n ${n} rel err ${err / scale}`);
    }
  }
});

test('Hines solve leaves A untouched and has a tiny residual', () => {
  const { A, b } = randomSystem(5000, 7, false);
  const diag0 = Float64Array.from(A.diag);
  const x = A.solve(Float64Array.from(b));
  assert.deepEqual(A.diag, diag0);
  const r = A.multiply(x, new Float64Array(A.n));
  let worst = 0;
  for (let i = 0; i < A.n; i++) worst = Math.max(worst, Math.abs(r[i]! - b[i]!));
  assert.ok(worst < 1e-12, `residual ${worst}`);
});

test('TreeMatrix rejects nodes that are not in Hines order', () => {
  assert.throws(() => new TreeMatrix(Int32Array.from([-1, 2, 0])), /Hines order/);
  assert.throws(() => new TreeMatrix(Int32Array.from([0, 0])), /root/);
});

test('leaves-first elimination creates no fill; other orders do', () => {
  const m = randomNodeTree(300, 3);
  const o = orders(m.n, rng(11));
  assert.equal(countFill(m.parent, o.hines), 0);
  assert.ok(countFill(m.parent, o.random) > 0);
  assert.ok(countFill(m.parent, o.rootFirst) > 0);
});

test('eliminating the hub of a star first fills the whole matrix', () => {
  const n = 40;
  const parent = new Int32Array(n).fill(0);
  parent[0] = -1;
  const rootFirst = Int32Array.from({ length: n }, (_, i) => i);
  assert.equal(countFill(parent, rootFirst), ((n - 1) * (n - 2)) / 2);
});

test('buildTree emits Hines order and the expected node counts', () => {
  const secs = randomSections(12, 5, 4);
  const withChildren = new Set(secs.slice(1).map((s) => s.parent)).size;
  const direct = buildTree(secs, { branch: 'direct' }).morph;
  const junction = buildTree(secs, { branch: 'junction', rootNode: true }).morph;
  assert.equal(direct.n, 12 * 5);
  assert.equal(junction.n, 12 * 5 + withChildren + 1);
  for (const m of [direct, junction]) {
    assert.equal(m.parent[0], -1);
    for (let i = 1; i < m.n; i++) assert.ok(m.parent[i]! < i && m.parent[i]! >= 0);
    new TreeMatrix(m.parent); // must not throw
  }
  const zeroArea = [...junction.length].filter((L) => L === 0).length;
  assert.equal(zeroArea, withChildren + 1);
});
