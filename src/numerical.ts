export interface PrimePoint {
  x: number;
  piX: number;
}

export interface SieveResult {
  counts: Record<number, number>;
  elapsedMs: number;
  bytesAllocated: number;
}

export type Matrix2x2 = [[number, number], [number, number]];
export type Vector2 = [number, number];

export interface NormalSystem2x2 {
  matrix: Matrix2x2;
  rhs: Vector2;
  sums: {
    n: number;
    sumU: number;
    sumY: number;
    sumU2: number;
    sumUY: number;
  };
}

export interface Gauss2x2Solution {
  params: Vector2;
  multiplier: number;
  upperMatrix: Matrix2x2;
  upperRhs: Vector2;
}

export interface LinearFitResult {
  a: number;
  b: number;
  normalSystem: NormalSystem2x2;
  gauss: Gauss2x2Solution;
  residualQ: number;
  predictions: Array<{
    x: number;
    u: number;
    realPiX: number;
    predictedPiX: number;
    deviation: number;
  }>;
  predict: (x: number) => number;
}

export interface LegendreFitResult {
  c0: number;
  c1: number;
  normalSystem: NormalSystem2x2;
  gauss: Gauss2x2Solution;
  residualQOriginalScale: number;
  residualQTransformedScale: number;
  predictions: Array<{
    x: number;
    X: number;
    wReal: number;
    wPredicted: number;
    realPiX: number;
    predictedPiX: number;
    deviationOriginal: number;
  }>;
  predict: (x: number) => number;
}

/**
 * Aula 01: Arredonda um número em ponto flutuante para `sigDigits` algarismos significativos.
 * Se `sigDigits` for omitido, mantém a precisão total de 64 bits (IEEE 754) do TypeScript.
 */
export function roundToSignificantDigits(value: number, sigDigits?: number): number {
  if (sigDigits === undefined || value === 0 || !Number.isFinite(value)) {
    return value;
  }
  return Number(value.toPrecision(sigDigits));
}

/**
 * Etapa 1 (E1): Crivo de Eratóstenes até `limit` usando `Uint8Array` (1 byte por número).
 */
export function sievePrimeCounts(limit: number, checkpoints: number[]): SieveResult {
  const start = performance.now();
  const isComposite = new Uint8Array(limit + 1);
  isComposite[0] = 1;
  isComposite[1] = 1;

  const sqrtLimit = Math.floor(Math.sqrt(limit));
  for (let p = 2; p <= sqrtLimit; p++) {
    if (isComposite[p] === 0) {
      for (let multiple = p * p; multiple <= limit; multiple += p) {
        isComposite[multiple] = 1;
      }
    }
  }

  const sortedCheckpoints = [...checkpoints].sort((a, b) => a - b);
  const counts: Record<number, number> = {};
  let runningCount = 0;
  let cpIndex = 0;

  for (let i = 2; i <= limit && cpIndex < sortedCheckpoints.length; i++) {
    if (isComposite[i] === 0) {
      runningCount++;
    }
    while (cpIndex < sortedCheckpoints.length && i === sortedCheckpoints[cpIndex]) {
      counts[sortedCheckpoints[cpIndex]] = runningCount;
      cpIndex++;
    }
  }

  const elapsedMs = performance.now() - start;
  return {
    counts,
    elapsedMs,
    bytesAllocated: isComposite.byteLength,
  };
}

/**
 * Aula 04: Monta o Sistema Normal 2x2 para o modelo linear `y ≈ param0 + param1 * u`.
 *
 * {   n * param0    + (Σ u_i) * param1 = Σ y_i
 * { (Σ u_i)*param0  + (Σ u_i²)* param1 = Σ (u_i * y_i)
 */
export function buildNormalSystem(
  u: number[],
  y: number[],
  sigDigits?: number,
): NormalSystem2x2 {
  const r = (val: number) => roundToSignificantDigits(val, sigDigits);
  const n = u.length;
  let sumU = 0;
  let sumY = 0;
  let sumU2 = 0;
  let sumUY = 0;

  for (let i = 0; i < n; i++) {
    const ui = r(u[i]);
    const yi = r(y[i]);
    sumU = r(sumU + ui);
    sumY = r(sumY + yi);
    sumU2 = r(sumU2 + r(ui * ui));
    sumUY = r(sumUY + r(ui * yi));
  }

  return {
    matrix: [
      [n, sumU],
      [sumU, sumU2],
    ],
    rhs: [sumY, sumUY],
    sums: { n, sumU, sumY, sumU2, sumUY },
  };
}

/**
 * Aula 03: Resolve um sistema linear 2x2 por Eliminação de Gauss explícita:
 * L2 <- L2 - m * L1, onde m = A[1][0] / A[0][0].
 */
export function solveGauss2x2(
  matrix: Matrix2x2,
  rhs: Vector2,
  sigDigits?: number,
): Gauss2x2Solution {
  const r = (val: number) => roundToSignificantDigits(val, sigDigits);
  const [[a00, a01], [a10, a11]] = matrix;
  const [b0, b1] = rhs;

  const multiplier = r(a10 / a00);
  const newA11 = r(a11 - r(multiplier * a01));
  const newB1 = r(b1 - r(multiplier * b0));

  const param1 = r(newB1 / newA11);
  const numerator0 = r(b0 - r(a01 * param1));
  const param0 = r(numerator0 / a00);

  return {
    params: [param0, param1],
    multiplier,
    upperMatrix: [
      [a00, a01],
      [0, newA11],
    ],
    upperRhs: [b0, newB1],
  };
}

/**
 * Etapa 2 (Ponte 1 / E2): Ajuste linear nos parâmetros `π(x) ≈ a + b * (x / log10(x))`.
 */
export function fitLinearModel(
  points: PrimePoint[],
  sigDigits?: number,
): LinearFitResult {
  const r = (val: number) => roundToSignificantDigits(val, sigDigits);
  const u = points.map((p) => r(p.x / Math.log10(p.x)));
  const y = points.map((p) => p.piX);

  const normalSystem = buildNormalSystem(u, y, sigDigits);
  const gauss = solveGauss2x2(normalSystem.matrix, normalSystem.rhs, sigDigits);
  const [a, b] = gauss.params;

  const predict = (x: number): number => a + b * (x / Math.log10(x));

  let residualQ = 0;
  const predictions = points.map((p, idx) => {
    const predictedPiX = predict(p.x);
    const deviation = predictedPiX - p.piX;
    residualQ += deviation * deviation;
    return {
      x: p.x,
      u: u[idx],
      realPiX: p.piX,
      predictedPiX,
      deviation,
    };
  });

  return {
    a,
    b,
    normalSystem,
    gauss,
    residualQ,
    predictions,
    predict,
  };
}

/**
 * Etapa 3 (Ponte 2 / E3): Forma de Legendre `π(x) ≈ x / (c0 + c1 * ln(x))`
 * resolvida por Linearização: `w = x / π(x) = c0 + c1 * X`, com `X = ln(x)`.
 */
export function fitLegendreModel(
  points: PrimePoint[],
  sigDigits?: number,
): LegendreFitResult {
  const r = (val: number) => roundToSignificantDigits(val, sigDigits);
  const X = points.map((p) => r(Math.log(p.x)));
  const w = points.map((p) => r(p.x / p.piX));

  const normalSystem = buildNormalSystem(X, w, sigDigits);
  const gauss = solveGauss2x2(normalSystem.matrix, normalSystem.rhs, sigDigits);
  const [c0, c1] = gauss.params;

  const predict = (x: number): number => x / (c0 + c1 * Math.log(x));

  let residualQOriginalScale = 0;
  let residualQTransformedScale = 0;

  const predictions = points.map((p, idx) => {
    const wPredicted = c0 + c1 * X[idx];
    const wDiff = wPredicted - w[idx];
    residualQTransformedScale += wDiff * wDiff;

    const predictedPiX = predict(p.x);
    const deviationOriginal = predictedPiX - p.piX;
    residualQOriginalScale += deviationOriginal * deviationOriginal;

    return {
      x: p.x,
      X: X[idx],
      wReal: w[idx],
      wPredicted,
      realPiX: p.piX,
      predictedPiX,
      deviationOriginal,
    };
  });

  return {
    c0,
    c1,
    normalSystem,
    gauss,
    residualQOriginalScale,
    residualQTransformedScale,
    predictions,
    predict,
  };
}

/**
 * Fórmula histórica de Legendre (1808): `c0 = -1.08366` e `c1 = 1`.
 */
export function predictLegendre1808(x: number): number {
  return x / (Math.log(x) - 1.08366);
}

export interface BenchmarkStats {
  warmups: number;
  runs: number;
  minMs: number;
  medianMs: number;
  meanMs: number;
}

/**
 * Executa um protocolo de benchmark pareado:
 * 1. Roda `warmups` vezes descartando os resultados (aquecimento do JIT / V8).
 * 2. Roda `runs` vezes medindo o tempo em milissegundos.
 * 3. Ordena as amostras e retorna Mínimo, Mediana e Média.
 */
export function runPairedBenchmark(
  fn: () => void,
  warmups = 5,
  runs = 30,
): BenchmarkStats {
  for (let i = 0; i < warmups; i++) {
    fn();
  }

  const samples: number[] = [];
  for (let i = 0; i < runs; i++) {
    const t0 = performance.now();
    fn();
    samples.push(performance.now() - t0);
  }

  samples.sort((a, b) => a - b);
  const minMs = samples[0];
  const mid = Math.floor(samples.length / 2);
  const medianMs =
    samples.length % 2 === 0 ? (samples[mid - 1] + samples[mid]) / 2 : samples[mid];
  const meanMs = samples.reduce((acc, v) => acc + v, 0) / samples.length;

  return {
    warmups,
    runs,
    minMs,
    medianMs,
    meanMs,
  };
}

/**
 * Etapa 4 (E4 & E5): Análise de chave RSA de 1024 bits (`x = 2^1024`) e limites do IEEE 754.
 */
export function analyzeRsa1024() {
  const directPow1024 = 2 ** 1024;
  const directOverflow = !Number.isFinite(directPow1024);

  // Solução via propriedade de logaritmo: ln(2^1024) = 1024 * ln(2)
  const ln2Pow1024 = 1024 * Math.LN2;
  const primeDensity = 1 / ln2Pow1024; // ~ 1 em cada 710 números
  const avgOddDraws = ln2Pow1024 / 2; // testando apenas ímpares: ~355 tentativas

  // Quantidade aproximada de primos de 1024 bits (entre 2^1023 e 2^1024):
  // N_primos ≈ 2^1023 / ln(2^1024). Calculamos em base 10 para não dar overflow:
  const log10TotalPrimes1024Bits = 1023 * Math.log10(2) - Math.log10(ln2Pow1024);
  const exponent10 = Math.floor(log10TotalPrimes1024Bits);
  const mantissa10 = 10 ** (log10TotalPrimes1024Bits - exponent10);

  return {
    directPow1024,
    directOverflow,
    ln2Pow1024,
    primeDensity,
    oneInEvery: Math.round(1 / primeDensity),
    avgOddDraws,
    scientificNotation1024BitPrimes: `${mantissa10.toFixed(3)} × 10^${exponent10}`,
    mantissa10,
    exponent10,
  };
}
