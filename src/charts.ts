import { mkdirSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { createCanvas } from '@napi-rs/canvas';
import { Chart, ChartConfiguration } from 'chart.js/auto';
import {
  fitLegendreModel,
  fitLinearModel,
  predictLegendre1808,
  PrimePoint,
  sievePrimeCounts,
} from './numerical.js';

const whiteBackgroundPlugin = {
  id: 'whiteBackground',
  beforeDraw: (chart: Chart) => {
    const { ctx, width, height } = chart;
    ctx.save();
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, width, height);
    ctx.restore();
  },
};

function renderChartToPng(config: ChartConfiguration, outputPath: string, width = 1200, height = 700) {
  const canvas = createCanvas(width, height);
  const ctx = canvas.getContext('2d');
  // Chart.js aceita o contexto 2D do @napi-rs/canvas diretamente
  const chart = new Chart(ctx as unknown as CanvasRenderingContext2D, {
    ...config,
    plugins: [whiteBackgroundPlugin, ...(config.plugins ?? [])],
    options: {
      ...config.options,
      animation: false,
      responsive: false,
    },
  });
  chart.update();
  const buffer = canvas.toBuffer('image/png');
  writeFileSync(outputPath, buffer);
  chart.destroy();
}

function generateCharts() {
  const outDir = join(process.cwd(), 'outputs');
  mkdirSync(outDir, { recursive: true });

  const checkpoints = [100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000];
  const sieve = sievePrimeCounts(10_000_000, checkpoints);

  const trainPoints: PrimePoint[] = [100, 1_000, 10_000, 100_000, 1_000_000].map((x) => ({
    x,
    piX: sieve.counts[x],
  }));

  const linearTrain = fitLinearModel(trainPoints);
  const legendreTrain = fitLegendreModel(trainPoints);

  // Gráfico 1: Erro Relativo (%) de cada aproximação em relação ao Crivo Real (10^2 a 10^7)
  const labels = ['10²', '10³', '10⁴', '10⁵', '10⁶ (treino)', '10⁷ (teste)'];
  const relErr = (pred: number, real: number) => Number((((pred - real) / real) * 100).toFixed(2));

  const errAsymptotic = checkpoints.map((x) => relErr(x / Math.log(x), sieve.counts[x]));
  const errLinear = checkpoints.map((x) => relErr(linearTrain.predict(x), sieve.counts[x]));
  const errLegendreFit = checkpoints.map((x) => relErr(legendreTrain.predict(x), sieve.counts[x]));
  const errLegendre1808 = checkpoints.map((x) =>
    relErr(predictLegendre1808(x), sieve.counts[x]),
  );

  const chart1Path = join(outDir, 'comparacao-pi-x.png');
  renderChartToPng(
    {
      type: 'line',
      data: {
        labels,
        datasets: [
          {
            label: 'Legendre Linearizado (Treino 10²..10⁶)',
            data: errLegendreFit,
            borderColor: '#1d4ed8',
            backgroundColor: '#1d4ed8',
            borderWidth: 3,
            pointRadius: 5,
          },
          {
            label: 'Mínimos Quadrados Linear (a + b·x/log₁₀x)',
            data: errLinear,
            borderColor: '#d97706',
            backgroundColor: '#d97706',
            borderWidth: 3,
            pointRadius: 5,
          },
          {
            label: 'Legendre Histórico (1808: c₀ = -1.08366)',
            data: errLegendre1808,
            borderColor: '#059669',
            backgroundColor: '#059669',
            borderWidth: 2,
            borderDash: [6, 4],
            pointRadius: 4,
          },
          {
            label: 'Assintótico Puro (x / ln x)',
            data: errAsymptotic,
            borderColor: '#dc2626',
            backgroundColor: '#dc2626',
            borderWidth: 2,
            borderDash: [3, 3],
            pointRadius: 4,
          },
        ],
      },
      options: {
        plugins: {
          title: {
            display: true,
            text: 'Erro Relativo (%) na Escala Original de π(x) — De 10² até 10⁷ (Teste)',
            font: { size: 20, weight: 'bold' },
          },
          legend: {
            labels: { font: { size: 14 } },
          },
        },
        scales: {
          x: {
            title: { display: true, text: 'x (escala log₁₀)', font: { size: 15, weight: 'bold' } },
          },
          y: {
            title: {
              display: true,
              text: 'Erro Relativo vs. Crivo Real (%)',
              font: { size: 15, weight: 'bold' },
            },
          },
        },
      },
    },
    chart1Path,
  );

  // Gráfico 2: Reta Linearizada de Legendre no plano (X = ln x, w = x / π(x))
  const allPoints: PrimePoint[] = checkpoints.map((x) => ({ x, piX: sieve.counts[x] }));
  const scatterReal = allPoints.map((p) => ({
    x: Number(Math.log(p.x).toFixed(3)),
    y: Number((p.x / p.piX).toFixed(4)),
  }));
  const lineFitted = allPoints.map((p) => {
    const X = Math.log(p.x);
    return {
      x: Number(X.toFixed(3)),
      y: Number((legendreTrain.c0 + legendreTrain.c1 * X).toFixed(4)),
    };
  });

  const chart2Path = join(outDir, 'reta-linearizada-legendre.png');
  renderChartToPng(
    {
      type: 'scatter',
      data: {
        datasets: [
          {
            label: 'Dados Reais Transformados: w = x / π(x)',
            data: scatterReal,
            backgroundColor: '#dc2626',
            pointRadius: 7,
          },
          {
            type: 'line',
            label: `Reta Ajustada (Mínimos Quadrados): w = ${legendreTrain.c0.toFixed(3)} + ${legendreTrain.c1.toFixed(3)}·ln(x)`,
            data: lineFitted,
            borderColor: '#1e3a8a',
            backgroundColor: '#1e3a8a',
            borderWidth: 3,
            pointRadius: 0,
          },
        ],
      },
      options: {
        plugins: {
          title: {
            display: true,
            text: 'Linearização da Forma de Legendre no Plano (X = ln x, w = x / π(x))',
            font: { size: 20, weight: 'bold' },
          },
          legend: {
            labels: { font: { size: 14 } },
          },
        },
        scales: {
          x: {
            type: 'linear',
            title: { display: true, text: 'X = ln(x)', font: { size: 15, weight: 'bold' } },
          },
          y: {
            title: { display: true, text: 'w = x / π(x)', font: { size: 15, weight: 'bold' } },
          },
        },
      },
    },
    chart2Path,
  );

  console.log(`Gráficos gerados com sucesso em:\n  - ${chart1Path}\n  - ${chart2Path}`);
}

generateCharts();
