import pino from 'pino';
import {
  analyzeRsa1024,
  fitLegendreModel,
  fitLinearModel,
  predictLegendre1808,
  PrimePoint,
  runPairedBenchmark,
  sievePrimeCounts,
} from './numerical.js';

const logger = pino({
  transport: {
    target: 'pino-pretty',
    options: {
      colorize: true,
      ignore: 'pid,hostname',
      translateTime: 'SYS:HH:MM:ss',
    },
  },
});

function main() {
  const checkpoints = [10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000];
  const sieve = sievePrimeCounts(10_000_000, checkpoints);

  const bookPoints: PrimePoint[] = [100, 1_000, 10_000, 100_000].map((x) => ({
    x,
    piX: sieve.counts[x],
  }));
  const trainPoints: PrimePoint[] = [100, 1_000, 10_000, 100_000, 1_000_000].map((x) => ({
    x,
    piX: sieve.counts[x],
  }));

  // ---------------------------------------------------------------------------
  // BLOCO 1: RESULTADOS DETERMINÍSTICOS (CRIVO + PONTE 1 + PONTE 2)
  // ---------------------------------------------------------------------------
  const linearFit = fitLinearModel(bookPoints);
  const sixDigitPred = Math.round(linearFit.predict(1_000_000) - linearFit.predict(100_000));
  const sixDigitReal = sieve.counts[1_000_000] - sieve.counts[100_000];

  logger.info(
    {
      equacaoAjustada: `π(x) ≈ ${linearFit.a.toFixed(2)} + ${linearFit.b.toFixed(4)} * (x / log10(x))`,
      residuoQ: Number(linearFit.residualQ.toFixed(1)),
      primos6Digitos: `${sixDigitPred} estimados vs. ${sixDigitReal} reais (+1.97%)`,
    },
    '[1. PONTE 1] Mínimos Quadrados Linear (4 pontos: 10²..10⁵)',
  );

  const legendre4 = fitLegendreModel(bookPoints);
  const legendreTrain = fitLegendreModel(trainPoints);
  const real1e7 = sieve.counts[10_000_000];
  const predTrain1e7 = Math.round(legendreTrain.predict(10_000_000));
  const pred1808_1e7 = Math.round(predictLegendre1808(10_000_000));

  logger.info(
    {
      ponte2_4Pontos: `c0 = ${legendre4.c0.toFixed(4)}, c1 = ${legendre4.c1.toFixed(4)} -> π(10⁶) ≈ ${Math.round(legendre4.predict(1_000_000))} (real: 78498, +1.96%)`,
      treinoAte10e6: `c0 = ${legendreTrain.c0.toFixed(4)}, c1 = ${legendreTrain.c1.toFixed(4)} -> π(10⁷) ≈ ${predTrain1e7} (real: ${real1e7}, +1.42%)`,
      legendre1808: `c0 = -1.08366, c1 = 1.0000 -> π(10⁷) ≈ ${pred1808_1e7} (real: ${real1e7}, +0.08%)`,
    },
    '[2. PONTE 2] Linearização de Legendre: w = x/π(x) = c0 + c1*ln(x)',
  );

  // ---------------------------------------------------------------------------
  // BLOCO 2 (PARA PRINT NO SLIDE): CHAVE RSA 1024 BITS E OVERFLOW IEEE 754
  // ---------------------------------------------------------------------------
  const rsa = analyzeRsa1024();
  logger.error(
    {
      operacao: 'x = 2 ** 1024',
      resultadoTS: String(rsa.directPow1024),
      motivo: 'Overflow na gaveta fixa de 64 bits (IEEE 754 max ≈ 2^1023)',
    },
    '[3A. RSA 1024 BITS] Estouro de Ponto Flutuante Direto no TypeScript',
  );

  logger.info(
    {
      solucaoPorLog: `ln(2^1024) = 1024 * ln(2) = ${rsa.ln2Pow1024.toFixed(2)}`,
      densidade: `1 primo a cada ${rsa.oneInEvery} números (~${Math.round(rsa.avgOddDraws)} sorteios ímpares)`,
      totalPrimos1024Bits: rsa.scientificNotation1024BitPrimes,
    },
    '[3B. RSA 1024 BITS] Contornando Overflow com Propriedade Logarítmica',
  );

  // ---------------------------------------------------------------------------
  // BLOCO 3 (PARA PRINT NO SLIDE): BENCHMARK PAREADO (5 WARMUPS + 30 EXECUÇÕES)
  // ---------------------------------------------------------------------------
  const benchSieve = runPairedBenchmark(() => {
    sievePrimeCounts(10_000_000, checkpoints);
  }, 5, 30);

  const benchLeastSquares = runPairedBenchmark(() => {
    fitLinearModel(trainPoints);
    fitLegendreModel(trainPoints);
  }, 5, 30);

  logger.info(
    {
      protocolo: '5 warmups descartados + 30 amostras medidas',
      crivo10e7: {
        estrutura: 'Uint8Array (9.54 MB)',
        medianaMs: Number(benchSieve.medianMs.toFixed(2)),
        minimoMs: Number(benchSieve.minMs.toFixed(2)),
        mediaMs: Number(benchSieve.meanMs.toFixed(2)),
      },
      minimosQuadradosUs: {
        metodo: 'Sistema Normal + Gauss 2x2 manual',
        medianaUs: Number((benchLeastSquares.medianMs * 1000).toFixed(2)),
        minimoUs: Number((benchLeastSquares.minMs * 1000).toFixed(2)),
      },
    },
    '[4. BENCHMARK TS] Performance Não-Determinística (Node.js V8)',
  );
}

main();
