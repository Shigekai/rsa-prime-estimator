import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import {
  analyzeRsa1024,
  fitLegendreModel,
  fitLinearModel,
  PrimePoint,
  runPairedBenchmark,
  sievePrimeCounts,
} from '../src/numerical.js';

describe('Projeto Cálculo Numérico — Contagem de Primos π(x) e RSA', () => {
  const checkpoints = [10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000];
  const sieve = sievePrimeCounts(10_000_000, checkpoints);

  it('Etapa 1 (E1): Crivo de Eratóstenes deve bater exatamente com os valores reais até 10^7', () => {
    assert.equal(sieve.counts[10], 4);
    assert.equal(sieve.counts[100], 25);
    assert.equal(sieve.counts[1_000], 168);
    assert.equal(sieve.counts[10_000], 1229);
    assert.equal(sieve.counts[100_000], 9592);
    assert.equal(sieve.counts[1_000_000], 78498);
    assert.equal(sieve.counts[10_000_000], 664579);
  });

  it('Etapa 2 (Ponte 1 / E2): Mínimos Quadrados Linear com 4 pontos deve bater com o gabarito', () => {
    const points: PrimePoint[] = [100, 1_000, 10_000, 100_000].map((x) => ({
      x,
      piX: sieve.counts[x],
    }));
    const fit = fitLinearModel(points);

    // Matriz Normal: [[4, 22883.333...], [22883.333..., 406363611.11...]]
    assert.equal(fit.normalSystem.matrix[0][0], 4);
    assert.ok(Math.abs(fit.normalSystem.matrix[0][1] - 22883.3333) < 0.01);
    assert.ok(Math.abs(fit.normalSystem.matrix[1][1] - 406363611.11) < 1);
    assert.deepEqual(fit.normalSystem.rhs, [11014, 194969750]);

    // Parâmetros a ≈ 12.825 (12.83) e b ≈ 0.479069 (0.4791)
    assert.ok(Math.abs(fit.a - 12.825) < 0.01);
    assert.ok(Math.abs(fit.b - 0.479069) < 0.0001);
    assert.ok(Math.abs(fit.residualQ - 506.3) < 0.5);

    // Primos de 6 dígitos: ~70263
    const sixDigitPrimes = fit.predict(1_000_000) - fit.predict(100_000);
    assert.equal(Math.round(sixDigitPrimes), 70263);
  });

  it('Etapa 3 (Ponte 2 / E3): Linearização de Legendre deve bater com 4 pontos e com treino/teste', () => {
    const points4: PrimePoint[] = [100, 1_000, 10_000, 100_000].map((x) => ({
      x,
      piX: sieve.counts[x],
    }));
    const fit4 = fitLegendreModel(points4);

    assert.ok(Math.abs(fit4.c0 - -0.3825) < 0.001);
    assert.ok(Math.abs(fit4.c1 - 0.93201) < 0.001);
    assert.ok(Math.abs(fit4.predict(1_000_000) - 80040) < 5);

    // Treino 10^2..10^6 e Teste em 10^7
    const trainPoints: PrimePoint[] = [100, 1_000, 10_000, 100_000, 1_000_000].map((x) => ({
      x,
      piX: sieve.counts[x],
    }));
    const fitTrain = fitLegendreModel(trainPoints);

    assert.ok(Math.abs(fitTrain.c0 - -0.53) < 0.01);
    assert.ok(Math.abs(fitTrain.c1 - 0.953) < 0.01);
    assert.ok(Math.abs(fitTrain.predict(10_000_000) - 674031) < 50);
  });

  it('Etapa 4 (E4 & E5): Análise de chave RSA de 1024 bits deve detectar overflow direto e resolver por log', () => {
    const rsa = analyzeRsa1024();
    assert.equal(rsa.directOverflow, true);
    assert.equal(rsa.directPow1024, Infinity);
    assert.ok(Math.abs(rsa.ln2Pow1024 - 709.7827) < 0.01);
    assert.equal(Math.round(rsa.avgOddDraws), 355);
    assert.equal(rsa.exponent10, 305);
  });

  it('Etapa 5 (Benchmark): Protocolo pareado deve descartar warmups e retornar Mediana e Mínimo válidos', () => {
    let calls = 0;
    const stats = runPairedBenchmark(() => {
      calls++;
    }, 5, 30);
    assert.equal(calls, 35);
    assert.equal(stats.warmups, 5);
    assert.equal(stats.runs, 30);
    assert.ok(stats.minMs >= 0);
    assert.ok(stats.medianMs >= stats.minMs);
  });
});
