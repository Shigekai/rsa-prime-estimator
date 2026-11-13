# RSA Prime Estimator — Cálculo Numérico (TypeScript × Julia)

Estimativa da função de contagem de primos `π(x)` e da densidade de primos em chaves criptográficas **RSA de 1024 bits (`x = 2¹⁰²⁴`)** usando o **Método dos Mínimos Quadrados** (caso linear nos parâmetros e linearização da Forma de Legendre), comparando **TypeScript (Node.js)** e **Julia**.

## Estrutura

- `src/numerical.ts`: implementações matemáticas puras (Crivo com `Uint8Array`, montagem do Sistema Normal 2×2, Eliminação de Gauss com controle de dígitos significativos, ajustes linear/linearizado e análise do RSA).
- `src/index.ts`: execução das 4 etapas em TypeScript com logs estruturados via `pino` + `pino-pretty`.
- `julia/main.jl`: execução espelhada das 4 etapas em Julia (`BitVector`, operador nativo `S \ y` e `BigFloat(2)^1024` sem overflow).
- `src/charts.ts`: geração dos gráficos em `outputs/*.png` com `chart.js` + `@napi-rs/canvas`.
- `test/numerical.test.ts`: suíte de testes validando os gabaritos das Pontes 1 e 2 e do Crivo.

## Comandos

```bash
npm run ts        # Executa as 4 etapas em TypeScript (logs pino-pretty)
npm run julia     # Executa as 4 etapas em Julia (BitVector, S \ y, BigFloat)
npm run graficos  # Gera os gráficos PNG na pasta outputs/
npm test          # Roda tsc --noEmit e os testes automatizados
```
