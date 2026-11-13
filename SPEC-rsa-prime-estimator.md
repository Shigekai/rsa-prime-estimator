# SPEC — Estimando π(x) e Primos em Chaves RSA (TypeScript × Julia)

**Status:** CONCLUÍDA (2026-11-12) — Implementação, testes, gráficos, benchmark pareado e relatório entregues.

## Resumo

Implementação prática para o projeto de Cálculo Numérico (UNESP Rio Claro) que estima a função de contagem de números primos `π(x)` a partir de dados empíricos gerados pelo Crivo de Eratóstenes até `10⁷`, aplicando o **Método dos Mínimos Quadrados** (no caso linear nos parâmetros e no caso não-linear via linearização da Forma de Legendre) e estendendo a análise para chaves criptográficas **RSA de 1024 bits (`x = 2¹⁰²⁴`)**.

O projeto executa o mesmo fluxo de 4 etapas em **TypeScript (Node.js)** e em **Julia**, permitindo comparar lado a lado o comportamento numérico de uma linguagem de propósito geral/web (baseada exclusivamente em `number` IEEE 754 de 64 bits e implementação manual da Eliminação de Gauss) com uma linguagem de computação científica (suporte nativo a `BitVector`, operador matricial `\` e ponto flutuante de precisão arbitrária `BigFloat`).

---

## Requisitos

1. **R1 — Fluxo espelhado de 4 etapas em TypeScript (`src/index.ts`) e Julia (`julia/main.jl`):**
   - **Etapa 1 (Crivo de Eratóstenes até `10⁷`):** Gerar os valores exatos de `π(10^k)` para `k = 1..7`, medindo tempo de execução e consumo teórico do array de crivo (`Uint8Array` de 1 byte/item no TS vs. `BitVector` de 1 bit/item na Julia).
   - **Etapa 2 (Mínimos Quadrados — Caso Linear nos Parâmetros / Ponte 1):** Ajustar o modelo `π(x) ≈ a + b · (x / log₁₀x)` nos 4 pontos do Problema 7.3 do livro (`10², 10³, 10⁴, 10⁵`), exibindo a matriz 2×2 do Sistema Normal, os parâmetros `a` e `b`, o resíduo `Q`, e a estimativa da quantidade de primos de 6 dígitos (`π(10⁶) − π(10⁵)`).
   - **Etapa 3 (Mínimos Quadrados — Linearização de Legendre / Ponte 2 + Treino/Teste):** Linearizar `π(x) ≈ x / (c₀ + c₁ · ln x)` para `x / π(x) = c₀ + c₁ · ln x`:
     - Rodar com os 4 pontos (`10²..10⁵`) para validar contra a Ponte 2 e prever `10⁶`.
     - Rodar o experimento de Treino (`10²..10⁶`) e Teste fora da amostra (`10⁷`), comparando os erros relativos na escala original de `π(x)` contra o valor real do Crivo e contra as constantes históricas de Legendre (`c₀ = -1.08366, c₁ = 1`).
   - **Etapa 4 (RSA de 1024 bits e Limites de Ponto Flutuante — Aulas 01 e 03):**
     - Demonstrar o impacto de baixa precisão/cancelamento subtrativo na matriz 2×2 da Ponte 1 (entradas de `4` até `4 × 10⁸`) e o limite de 53 bits da mantissa (`2**53 + 1 === 2**53`).
     - Avaliar o caso de uma chave RSA de 1024 bits (`x = 2¹⁰²⁴`): demonstrar o **Overflow** (`Infinity`) na gaveta de 64 bits (`number` / `Float64`), resolver em TypeScript via propriedade logarítmica (`ln(2¹⁰²⁴) = 1024 · ln 2 ≈ 709.78`) e resolver diretamente em Julia via gaveta elástica (`BigFloat(2)^1024`), calculando a densidade de primos (`~1/710`), a média de sorteios ímpares (`~355 tentativas`) e o total de primos de 1024 bits (`~1.3 × 10³⁰⁵`).

2. **R2 — Módulo numérico puro e testável em TypeScript (`src/numerical.ts`):**
   - Funções puras e tipadas para o Crivo (`sievePrimeCounts`), montagem do Sistema Normal 2×2 (`buildNormalSystem`), resolução por Eliminação de Gauss 2×2 com opção de arredondamento para k dígitos significativos (`solveGauss2x2`), ajuste linear (`fitLinearModel`), ajuste linearizado de Legendre (`fitLegendreModel`) e cálculo do resíduo `Q`.

3. **R3 — Logs visuais prontos para prints de apresentação:**
   - Usar `pino` + `pino-pretty` no TypeScript com formatação limpa e níveis semânticos (`INFO` para resultados, `WARN`/`ERROR` ao demonstrar cancelamento subtrativo e overflow `Infinity`), facilitando a captura de screenshots para os slides de 10 minutos.

4. **R4 — Geração automatizada de gráficos (`src/charts.ts`):**
   - Usar `chart.js` + `@napi-rs/canvas` para gerar e salvar em `outputs/`:
     - `outputs/comparacao-pi-x.png`: curva real `π(x)` vs. modelos ajustados em escala logarítmica.
     - `outputs/reta-linearizada-legendre.png`: pontos `(ln x, x/π(x))` e a reta ajustada da Ponte 2.

5. **R5 — Interface enxuta de comandos no `package.json`:**
   - Apenas 4 scripts principais:
     - `npm run ts`: executa o fluxo completo de 4 etapas em TypeScript.
     - `npm run julia`: executa o fluxo completo de 4 etapas em Julia.
     - `npm run graficos`: gera os gráficos em `outputs/`.
     - `npm test`: executa `tsc --noEmit` e a suíte de testes automatizados (`node --import tsx --test`).

---

## Critérios de Aceite

- [x] **CA1:** `npm test` passa sem erros (`tsc --noEmit` limpo + testes unitários validando:
  - Crivo: `π(10^k)` para `k = 1..7` igual a `[4, 25, 168, 1229, 9592, 78498, 664579]`.
  - Ponte 1 (4 pontos): `a ≈ 12.825`, `b ≈ 0.479069`, resíduo `Q ≈ 506.3`, primos de 6 dígitos `≈ 70263`.
  - Ponte 2 (4 pontos): `c₀ ≈ -0.3825`, `c₁ ≈ 0.93201`, previsão `π(10⁶) ≈ 80040`.
  - Ponte 2 (treino `10²..10⁶`): `c₀ ≈ -0.530`, `c₁ ≈ 0.953`, previsão `π(10⁷) ≈ 674031`).
- [x] **CA2:** `npm run ts` executa as 4 etapas e exibe os logs formatados via `pino-pretty` sem falhas.
- [x] **CA3:** `npm run julia` executa as mesmas 4 etapas em Julia (`julia v1.13.3`) com resultados consistentes com o TypeScript e cálculo exato em `BigFloat` para `2¹⁰²⁴`.
- [x] **CA4:** `npm run graficos` cria os arquivos `outputs/comparacao-pi-x.png` e `outputs/reta-linearizada-legendre.png`.

---

## Restrições e Premissas

- **Ambiente local:** macOS (Apple Silicon), Node.js `v24.22.2`, `npm`, Julia `v1.13.3` (`/opt/homebrew/bin/julia`).
- **Sem frameworks web:** Sem NestJS/Express para não adicionar complexidade irrelevante à avaliação de Cálculo Numérico.
- **Escopo da matéria:** Foco exclusivo nos conteúdos já estudados nas Pontes 1 e 2 (Aula 04 — Mínimos Quadrados Linear e Linearização) somados à Aula 01 (Ponto Flutuante) e Aula 03 (Eliminação de Gauss). Interpolação (Aula 05) e Integração (Aula 06) ficam fora do escopo do projeto.

---

## Abordagem Técnica (Plano Aprovado — Fase 2)

1. **Setup (`package.json`, `tsconfig.json`, `.gitignore`):**
   - Configurar TypeScript (`strict: true`, ESM) e instalar `pino`, `pino-pretty`, `chart.js`, `@napi-rs/canvas`, `typescript`, `tsx`, `@types/node`.
   - Definir os 4 comandos no `package.json`: `ts`, `julia`, `graficos`, `test`.
2. **Núcleo Matemático (`src/numerical.ts`):**
   - `sievePrimeCounts`: Crivo de Eratóstenes com `Uint8Array` até `10⁷`.
   - `roundToSignificantDigits`: arredondamento para `k` algarismos significativos (Aula 01).
   - `buildNormalSystem` e `solveGauss2x2`: montagem do Sistema Normal 2×2 e Eliminação de Gauss explícita com suporte opcional a precisão reduzida.
   - `fitLinearModel` (Ponte 1) e `fitLegendreModel` (Ponte 2), calculando predições e resíduo `Q` na escala original de `π(x)`.
   - `analyzeRsa1024`: demonstração do overflow `2 ** 1024 === Infinity` e cálculo via `1024 * Math.LN2`.
3. **Runner TypeScript (`src/index.ts`):**
   - Executar as 4 etapas com logs estruturados via `pino` + `pino-pretty` (`INFO` para resultados e `WARN`/`ERROR` para cancelamento subtrativo e overflow).
4. **Runner Julia (`julia/main.jl`):**
   - Executar as mesmas 4 etapas usando `LinearAlgebra` e `Printf` nativos (`BitVector` no crivo, `S \ y` nos ajustes e `BigFloat(2)^1024` sem overflow no RSA).
5. **Gráficos (`src/charts.ts`):**
   - Gerar `outputs/comparacao-pi-x.png` e `outputs/reta-linearizada-legendre.png` via `chart.js` + `@napi-rs/canvas`.
6. **Testes e Docs (`test/numerical.test.ts`, `README.md`):**
   - Validar todos os valores de referência com `node:test` + `tsc --noEmit`.

---

## Questões em Aberto

Nenhuma questão em aberto.

---

## Log de Decisões

- **2026-11-12:** Escopo enxuto definido: manter o tema de contagem de primos `π(x)` e RSA focando em Mínimos Quadrados (Ponte 1 linear + Ponte 2 linearização de Legendre) + comparativo numérico TypeScript × Julia (Aulas 01 e 03), movendo Interpolação e Integração para estudo posterior da P2.
- **2026-11-12:** Descartado o uso de NestJS em favor de scripts diretos em TypeScript (`tsx`) e Julia (`julia 1.13.3`), priorizando clareza matemática e gestão de tempo.
- **2026-11-12:** Definido uso de `pino` + `pino-pretty` no TypeScript para gerar logs estruturados e coloridos que servirão de prints nos slides da apresentação de 10 minutos.
- **2026-11-12:** Definido que `npm run ts` e `npm run julia` executarão exatamente o mesmo fluxo de 4 etapas (Crivo → Mínimos Quadrados Linear → Linearização de Legendre → RSA 1024 bits e limites de ponto flutuante).
- **2026-11-12:** Reduzida a lista de scripts do `package.json` para apenas 4 comandos essenciais: `ts`, `julia`, `graficos` e `test`.
- **2026-11-12:** Geração dos gráficos `.png` será feita via `chart.js` + `@napi-rs/canvas` em `src/charts.ts`.
- **2026-11-12:** Spec aprovada pelo usuário (fim da Fase 1) e plano de implementação aprovado pelo usuário (fim da Fase 2). Início da Fase 3.
- **2026-11-12:** Adicionado protocolo de benchmark pareado não-determinístico (5 rodadas de warmup descartadas + 30 amostras medidas reportando Mediana, Mínimo e Média) no TypeScript e na Julia, além da consolidação do relatório completo e tabelas determinísticas/não-determinísticas em `RELATORIO.md`.
- **2026-11-12:** Removido bloco de cancelamento subtrativo do log principal a pedido do usuário para manter o foco nas 4 etapas. Revisão humana aprovada e tarefa concluída (Fase 4).
