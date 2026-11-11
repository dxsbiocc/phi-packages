/**
 * Small, dependency-free statistics for biomarker associations.
 *
 * These are screening statistics over cell-line panels. They carry no multiple-testing
 * correction and treat models as independent, so callers report them as descriptive.
 */

const EPS = 3e-14;
const TINY = 1e-300;

const LANCZOS = [
  57.1562356658629235, -59.5979603554754912, 14.1360979747417471, -0.491913816097620199,
  0.339946499848118887e-4, 0.465236289270485756e-4, -0.983744753048795646e-4,
  0.158088703224912494e-3, -0.210264441724104883e-3, 0.217439618115212643e-3,
  -0.164318106536763890e-3, 0.844182239838527433e-4, -0.261908384015814087e-4,
  0.368991826595316234e-5,
] as const;

function logGamma(x: number): number {
  let y = x;
  const tmp = x + 5.24218750000000000;
  const base = (x + 0.5) * Math.log(tmp) - tmp;
  let series = 0.999999999999997092;
  for (const coefficient of LANCZOS) {
    y += 1;
    series += coefficient / y;
  }
  return base + Math.log((2.5066282746310005 * series) / x);
}

function betaContinuedFraction(a: number, b: number, x: number): number {
  const qab = a + b;
  const qap = a + 1;
  const qam = a - 1;
  let c = 1;
  let d = 1 - (qab * x) / qap;
  if (Math.abs(d) < TINY) d = TINY;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= 500; m += 1) {
    const m2 = 2 * m;
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < TINY) d = TINY;
    c = 1 + aa / c;
    if (Math.abs(c) < TINY) c = TINY;
    d = 1 / d;
    h *= d * c;
    aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < TINY) d = TINY;
    c = 1 + aa / c;
    if (Math.abs(c) < TINY) c = TINY;
    d = 1 / d;
    const delta = d * c;
    h *= delta;
    if (Math.abs(delta - 1) < EPS) break;
  }
  return h;
}

/** Regularized incomplete beta function I_x(a, b). */
export function regularizedBeta(x: number, a: number, b: number): number {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const front = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  return x < (a + 1) / (a + b + 2)
    ? (front * betaContinuedFraction(a, b, x)) / a
    : 1 - (front * betaContinuedFraction(b, a, 1 - x)) / b;
}

/** Two-sided p-value of a normal deviate z. */
export function normalTwoSidedP(z: number): number {
  // erfc(|z|/sqrt 2) = Q(1/2, z^2/2), the regularized upper incomplete gamma function.
  const x = (z * z) / 2;
  if (x === 0) return 1;
  const a = 0.5;
  if (x < a + 1) {
    let sum = 1 / a;
    let term = sum;
    for (let n = 1; n < 500; n += 1) {
      term *= x / (a + n);
      sum += term;
      if (Math.abs(term) < Math.abs(sum) * EPS) break;
    }
    return Math.min(1, Math.max(0, 1 - sum * Math.exp(-x + a * Math.log(x) - logGamma(a))));
  }
  let b = x + 1 - a;
  let c = 1 / TINY;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i < 500; i += 1) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < TINY) d = TINY;
    c = b + an / c;
    if (Math.abs(c) < TINY) c = TINY;
    d = 1 / d;
    const delta = d * c;
    h *= delta;
    if (Math.abs(delta - 1) < EPS) break;
  }
  return Math.min(1, Math.max(0, Math.exp(-x + a * Math.log(x) - logGamma(a)) * h));
}

/** Average ranks (1-based); tied values share the mean of the ranks they span. */
export function averageRanks(values: readonly number[]): number[] {
  const order = values.map((value, index) => ({ value, index })).sort((left, right) => left.value - right.value);
  const ranks = new Array<number>(values.length);
  let start = 0;
  while (start < order.length) {
    let end = start;
    while (end + 1 < order.length && order[end + 1]!.value === order[start]!.value) end += 1;
    const rank = (start + end) / 2 + 1;
    for (let k = start; k <= end; k += 1) ranks[order[k]!.index] = rank;
    start = end + 1;
  }
  return ranks;
}

export function median(values: readonly number[]): number | null {
  if (!values.length) return null;
  const sorted = [...values].sort((left, right) => left - right);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle]! : (sorted[middle - 1]! + sorted[middle]!) / 2;
}

export function pearson(x: readonly number[], y: readonly number[]): number | null {
  if (x.length !== y.length || x.length < 2) return null;
  const mx = x.reduce((sum, value) => sum + value, 0) / x.length;
  const my = y.reduce((sum, value) => sum + value, 0) / y.length;
  let sxy = 0;
  let sxx = 0;
  let syy = 0;
  for (let i = 0; i < x.length; i += 1) {
    const dx = x[i]! - mx;
    const dy = y[i]! - my;
    sxy += dx * dy;
    sxx += dx * dx;
    syy += dy * dy;
  }
  if (sxx === 0 || syy === 0) return null;
  return sxy / Math.sqrt(sxx * syy);
}

export interface Correlation {
  n: number;
  rho: number;
  p: number;
}

/** Spearman rank correlation with a Student-t two-sided p-value (needs n >= 3 and varying data). */
export function spearman(x: readonly number[], y: readonly number[]): Correlation | null {
  const n = x.length;
  if (n !== y.length || n < 3) return null;
  const rho = pearson(averageRanks(x), averageRanks(y));
  if (rho == null) return null;
  const clamped = Math.max(-1, Math.min(1, rho));
  if (Math.abs(clamped) === 1) return { n, rho: clamped, p: 0 };
  const df = n - 2;
  const t = clamped * Math.sqrt(df / (1 - clamped * clamped));
  return { n, rho: clamped, p: regularizedBeta(df / (df + t * t), df / 2, 0.5) };
}

export interface RankTest {
  n1: number;
  n2: number;
  u: number;
  p: number;
  /** Rank-biserial correlation: positive means the first group tends to have larger values. */
  effect: number;
}

/** Mann-Whitney U with tie correction and continuity correction, two-sided normal approximation. */
export function mannWhitney(first: readonly number[], second: readonly number[]): RankTest | null {
  const n1 = first.length;
  const n2 = second.length;
  if (n1 === 0 || n2 === 0) return null;
  const pooled = [...first, ...second];
  const ranks = averageRanks(pooled);
  let rankSum = 0;
  for (let i = 0; i < n1; i += 1) rankSum += ranks[i]!;
  const u = rankSum - (n1 * (n1 + 1)) / 2;
  const n = n1 + n2;
  const counts = new Map<number, number>();
  for (const value of pooled) counts.set(value, (counts.get(value) ?? 0) + 1);
  let tieTerm = 0;
  for (const count of counts.values()) tieTerm += count ** 3 - count;
  const variance = ((n1 * n2) / 12) * (n + 1 - (n > 1 ? tieTerm / (n * (n - 1)) : 0));
  const effect = (2 * u) / (n1 * n2) - 1;
  if (variance <= 0) return { n1, n2, u, p: 1, effect };
  const mean = (n1 * n2) / 2;
  const deviation = Math.abs(u - mean);
  const z = Math.max(0, deviation - 0.5) / Math.sqrt(variance);
  return { n1, n2, u, p: normalTwoSidedP(z), effect };
}
