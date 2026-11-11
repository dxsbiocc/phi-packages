import assert from "node:assert/strict";
import test from "node:test";
import { averageRanks, mannWhitney, median, normalTwoSidedP, pearson, spearman } from "../stats.js";

// Reference values were computed with scipy 1.14 (scipy.stats.spearmanr, mannwhitneyu, norm.sf).
function close(actual: number, expected: number, relative = 1e-6): void {
  assert.ok(Math.abs(actual - expected) <= Math.abs(expected) * relative + 1e-12, `${actual} is not close to ${expected}`);
}

test("average ranks share the mean rank between ties", () => {
  assert.deepEqual(averageRanks([10, 20, 20, 30]), [1, 2.5, 2.5, 4]);
  assert.deepEqual(averageRanks([3, 1, 2]), [3, 1, 2]);
});

test("normal two-sided p-values match scipy", () => {
  close(normalTwoSidedP(1.96), 0.04999579029644087);
  close(normalTwoSidedP(0.3), 0.7641771556220948);
  close(normalTwoSidedP(6), 1.9731752900753887e-9, 1e-4);
  assert.equal(normalTwoSidedP(0), 1);
});

test("spearman correlation and p-value match scipy, including ties", () => {
  const rho = spearman([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [2, 1, 4, 3, 7, 8, 5, 10, 9, 12]);
  assert.ok(rho);
  close(rho.rho, 0.9272727272727272);
  close(rho.p, 0.00011203450639397582, 1e-5);
  const tied = spearman([1, 2, 2, 3, 4, 4, 4, 5], [5, 3, 4, 2, 1, 2, 2, 0]);
  assert.ok(tied);
  close(tied.rho, -0.9434148851749269);
  close(tied.p, 0.00043394122476300305, 1e-5);
});

test("spearman refuses too few points and constant data", () => {
  assert.equal(spearman([1, 2], [1, 2]), null);
  assert.equal(spearman([1, 1, 1, 1], [1, 2, 3, 4]), null);
});

test("mann-whitney matches scipy with ties and reports the direction of the effect", () => {
  const test1 = mannWhitney([0.2, 1.0, 1.1, 1.5, 2.0, 2.0], [2.5, 3.0, 3.5, 2.0, 4.0]);
  assert.ok(test1);
  assert.equal(test1.u, 1);
  close(test1.p, 0.012865581294556065, 1e-5);
  assert.ok(test1.effect < 0, "the first group has the smaller values");
  close(test1.effect, 2 / 30 - 1);
  const swapped = mannWhitney([2.5, 3.0, 3.5, 2.0, 4.0], [0.2, 1.0, 1.1, 1.5, 2.0, 2.0]);
  assert.ok(swapped);
  close(swapped.p, test1.p);
  close(swapped.effect, -test1.effect);
  assert.equal(mannWhitney([], [1]), null);
});

test("pearson and median handle the basic cases", () => {
  close(pearson([1, 2, 3], [2, 4, 6])!, 1);
  assert.equal(pearson([1, 1, 1], [1, 2, 3]), null);
  assert.equal(median([3, 1, 2]), 2);
  assert.equal(median([4, 1, 2, 3]), 2.5);
  assert.equal(median([]), null);
});
