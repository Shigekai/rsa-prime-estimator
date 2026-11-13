# Quantos primos cabem numa chave RSA? Estimando π(x) por Mínimos Quadrados e Linearização (TypeScript × Julia)

**Autor:** Lucas · Ciência da Computação (UNESP Rio Claro)  
**Disciplina:** Cálculo Numérico — Prof.ª Helena Morais  
**Repositório:** `https://github.com/Shigekai/rsa-prime-estimator`

---

## Guia Rápido para Montar os Slides (Determinístico vs. Não-Determinístico)

### A. O que colocar como TABELA NATIVA nos slides (Dados Determinísticos)
Não use print de terminal para matemática determinística — copie as tabelas das seções abaixo direto para os slides:
1. **Tabela 1 (Crivo até `10⁷`):** valores exatos de `π(x)` (`4, 25, 168, 1229, 9592, 78498, 664579`).
2. **Tabela 2 (Ponte 1 — Mínimos Quadrados Linear):** matriz 2×2 do Sistema Normal, `a = 12.83`, `b = 0.4791`, `Q = 506.3`, e a previsão de primos de 6 dígitos (`70 263` vs. `68 906` reais, erro `+1.97%`).
3. **Tabela 3 (Ponte 2 — Linearização de Legendre e Campeonato em `10⁷`):** comparação de erro relativo na escala original de `π(x)`.
4. **Gráficos `.png` gerados em `outputs/`:**
   - `outputs/reta-linearizada-legendre.png` (mostra os pontos `(ln x, x/π(x))` virando uma reta perfeita).
   - `outputs/comparacao-pi-x.png` (mostra a queda do erro relativo até `10⁷`).

### B. O que colocar como PRINT DE TERMINAL nos slides (Só 2 prints!)
Reserve prints do terminal apenas para mostrar **comportamento de execução das linguagens**:
- **Print 1 (Slide 6 — Ponto Flutuante e Overflow no RSA):** Os blocos `[3A. RSA 1024 BITS]` e `[3B. RSA 1024 BITS]` do `npm run ts` lado a lado com o bloco `[3. RSA 1024 BITS]` do `npm run julia`. O badge vermelho (`ERROR: Infinity` no `2 ** 1024`) seguido da solução por logaritmo prova visualmente os limites do IEEE 754 de 64 bits frente ao `BigFloat` da Julia.
- **Print 2 (Slide 6/7 — Benchmark Não-Determinístico):** Os blocos `[4. BENCHMARK TS]` e `[4. BENCHMARK JULIA]` mostrando o protocolo de 5 warmups + 30 amostras (Mediana e Mínimo).

---
---

# Relatório do Projeto (4 Seções Obrigatórias)

## 1. Introdução

A função de contagem de números primos, denotada por $\pi(x)$, fornece a quantidade de números primos menores ou iguais a um valor real positivo $x$. Embora sua definição seja elementar, não existe uma fórmula algébrica fechada eficiente para calcular $\pi(x)$ diretamente para valores grandes de $x$. Na prática computacional, determinar $\pi(x)$ de forma exata exige peneirar todos os inteiros até $x$ (por exemplo, via Crivo de Eratóstenes, com custo de tempo $O(x \ln \ln x)$ e memória $O(x)$), o que se torna inviável em escalas criptográficas.

No final do século XVIII e início do século XIX, antes da existência de computadores eletrônicos, Carl Friedrich Gauss (~1792) e Adrien-Marie Legendre (1808) inspecionaram tabelas de números primos calculadas à mão e observaram que a densidade de primos em torno de $x$ decai de forma inversamente proporcional ao logaritmo natural $\ln(x)$. Legendre propôs a aproximação racional $\pi(x) \approx \frac{x}{c_0 + c_1 \ln(x)}$ (estimando empiricamente $c_1 = 1$ e $c_0 = -1.08366$), relação que mais tarde culminou na demonstração do Teorema dos Números Primos ($\pi(x) \sim \frac{x}{\ln x}$) em 1896.

Na Ciência da Computação moderna, essa estimativa assintótica é a base da geração de chaves no algoritmo criptográfico **RSA**: para construir uma chave de **1024 bits** (números da ordem de $x = 2^{1024} \approx 1.79 \times 10^{308}$), o sistema sorteia números ímpares aleatórios e testa sua primalidade. Saber quantos primos existem nessa faixa determina quantas tentativas, em média, o algoritmo precisa realizar até encontrar um número primo.

**Objetivo:** Este trabalho aplica o **Método dos Mínimos Quadrados** (no caso linear nos parâmetros e no caso não-linear via linearização da Forma de Legendre) sobre dados reais gerados pelo Crivo de Eratóstenes até $x = 10^7$ para estimar a lei de formação de $\pi(x)$. Adicionalmente, implementamos o mesmo fluxo em **TypeScript (Node.js)** e em **Julia**, analisando tanto os resultados matemáticos determinísticos (incluindo cancelamento subtrativo e condicionamento do Sistema Normal) quanto os limites de ponto flutuante (Overflow do padrão IEEE 754 de 64 bits em $2^{1024}$ frente à precisão arbitrária `BigFloat`) e um benchmark pareado não-determinístico de tempo e memória.

---

## 2. Desenvolvimento

### 2.1 Geração dos Dados Empíricos: Crivo de Eratóstenes (Determinístico)

Para construir a tabela base de valores reais de $\pi(x)$, implementamos o Crivo de Eratóstenes até $x = 10^7$ nas duas linguagens. A contagem obtida (Tabela 1) é 100% determinística e serve como gabarito ("verdade terrestre") para avaliar os ajustes numéricos.

**Tabela 1 — Contagem exata de números primos $\pi(x)$ gerada pelo Crivo de Eratóstenes:**

| $x$ | $10^1$ | $10^2$ | $10^3$ | $10^4$ | $10^5$ | $10^6$ | $10^7$ |
|---|---|---|---|---|---|---|---|
| **$\pi(x)$ real** | 4 | 25 | 168 | 1 229 | 9 592 | 78 498 | 664 579 |

---

### 2.2 Mínimos Quadrados: Caso Linear nos Parâmetros (Ponte 1 / Problema 7.3)

Seguindo o Problema 7.3 de Franco (2006), buscamos ajustar aos $n = 4$ pontos iniciais ($x = 10^2, 10^3, 10^4, 10^5$) o modelo:

$$\pi(x) \approx a + b \cdot \frac{x}{\log_{10}(x)}$$

Como os parâmetros $a$ e $b$ possuem expoente $1$ no numerador, o modelo é **linear nos parâmetros**, com funções-base $u_0(x) = 1$ e $u_1(x) = u = \frac{x}{\log_{10}(x)}$. Minimizando o resíduo quadrático $Q(a, b) = \sum_{i=1}^{4} (a + b \cdot u_i - y_i)^2$ através de $\frac{\partial Q}{\partial a} = 0$ e $\frac{\partial Q}{\partial b} = 0$, obtemos o **Sistema Normal 2×2**:

$$\begin{bmatrix} n & \sum u_i \\ \sum u_i & \sum u_i^2 \end{bmatrix} \begin{bmatrix} a \\ b \end{bmatrix} = \begin{bmatrix} \sum y_i \\ \sum u_i y_i \end{bmatrix}$$

**Tabela 2 — Somatórios e Sistema Normal do ajuste linear nos parâmetros ($10^2$ a $10^5$):**

| Ponto $i$ | $x_i$ | $u_i = x_i / \log_{10}(x_i)$ | $y_i = \pi(x_i)$ | $u_i^2$ | $u_i \cdot y_i$ | $\hat{y}_i$ (Ajuste) | Desvio ($\hat{y}_i - y_i$) |
|---|---|---|---|---|---|---|---|
| 1 | $10^2$ | $50.0000$ | $25$ | $2\,500.0$ | $1\,250$ | $36.78$ | $+11.78$ |
| 2 | $10^3$ | $333.3333$ | $168$ | $111\,111.1$ | $56\,000$ | $172.51$ | $+4.51$ |
| 3 | $10^4$ | $2\,500.0000$ | $1\,229$ | $6\,250\,000.0$ | $3\,072\,500$ | $1\,210.50$ | $-18.50$ |
| 4 | $10^5$ | $20\,000.0000$ | $9\,592$ | $400\,000\,000.0$ | $191\,840\,000$ | $9\,594.21$ | $+2.21$ |
| **Soma ($\Sigma$)** | — | **$22\,883.3333$** | **$11\,014$** | **$406\,363\,611.1$** | **$194\,969\,750$** | — | **$Q = 506.32$** |

Resolvendo o sistema $\begin{bmatrix} 4 & 22\,883.333 \\ 22\,883.333 & 406\,363\,611.1 \end{bmatrix} \begin{bmatrix} a \\ b \end{bmatrix} = \begin{bmatrix} 11\,014 \\ 194\,969\,750 \end{bmatrix}$ por Eliminação de Gauss ($L_2 \leftarrow L_2 - 5\,720.8333 \cdot L_1$), tanto a implementação manual em TypeScript quanto o operador nativo `S \ y` da Julia convergem para:

- **$a = 12.8253 \approx 12.83$**
- **$b = 0.479069 \approx 0.4791$**
- **Resíduo quadrático:** $Q = 506.32$

**Estimativa de primos de 6 dígitos ($\pi(10^6) - \pi(10^5)$):**
Aplicando o modelo ajustado, obtemos $\pi(10^6) \approx 79\,857.7$ e $\pi(10^5) \approx 9\,594.2$, resultando em **$70\,263$ primos de 6 dígitos**. Comparado ao valor real do Crivo ($78\,498 - 9\,592 = 68\,906$), o erro relativo é de apenas **$+1.97\%$**. Além disso, convertendo a base do logaritmo ($b \cdot \frac{x}{\log_{10}x} = (b \ln 10) \cdot \frac{x}{\ln x}$), encontramos $b \ln 10 = 0.479069 \times 2.302585 = \mathbf{1.1031}$, redescobrindo empiricamente o coeficiente unitário do Teorema dos Números Primos.

---

### 2.3 Mínimos Quadrados: Linearização da Forma de Legendre (Ponte 2)

A Forma de Legendre propõe o modelo:

$$\pi(x) \approx \frac{x}{c_0 + c_1 \ln(x)}$$

Neste modelo, os parâmetros $c_0$ e $c_1$ encontram-se no **denominador** (expoente $-1$), o que impede a montagem direta de uma matriz linear. Para aplicar o mesmo método do Sistema Normal, realizamos o pré-processamento algébrico de **linearização**, invertendo a relação:

$$w = \frac{x}{\pi(x)} = c_0 + c_1 \ln(x)$$

Definindo as variáveis transformadas $X_i = \ln(x_i)$ e $w_i = \frac{x_i}{\pi(x_i)}$, os parâmetros voltam a ter expoente $+1$ no numerador ($w = c_0 + c_1 X$), permitindo aplicar o mesmo Sistema Normal 2×2.

**Tabela 3 — Dados transformados da Linearização de Legendre ($X_i = \ln x_i$, $w_i = x_i / \pi(x_i)$):**

| Ponto $i$ | $x_i$ | $\pi(x_i)$ real | $X_i = \ln(x_i)$ | $w_i = x_i / \pi(x_i)$ | $X_i^2$ | $X_i \cdot w_i$ |
|---|---|---|---|---|---|---|
| 1 | $10^2$ | $25$ | $4.60517$ | $4.00000$ | $21.2076$ | $18.4207$ |
| 2 | $10^3$ | $168$ | $6.90776$ | $5.95238$ | $47.7171$ | $41.1176$ |
| 3 | $10^4$ | $1\,229$ | $9.21034$ | $8.13670$ | $84.8304$ | $74.9417$ |
| 4 | $10^5$ | $9\,592$ | $11.51293$ | $10.42535$ | $132.5475$ | $120.0263$ |
| **Soma (4 pts)** | — | — | **$32.23620$** | **$28.51443$** | **$286.3025$** | **$254.5063$** |
| 5 (Treino) | $10^6$ | $78\,498$ | $13.81551$ | $12.73918$ | $190.8683$ | $176.0003$ |

Realizamos dois experimentos determinísticos com a forma linearizada:
1. **Ajuste com 4 pontos ($10^2..10^5$, Ponte 2):** O Sistema Normal fornece **$c_0 = -0.3825$** e **$c_1 = 0.93201$**. Devolvendo os parâmetros para a escala original $\pi(x) = \frac{x}{-0.3825 + 0.93201 \ln x}$, a extrapolação para $x = 10^6$ prevê **$80\,040$ primos** (contra $78\,498$ reais, erro de **$+1.96\%$**).
2. **Treino ($10^2..10^6$) e Teste fora da amostra ($10^7$):** Ajustando o modelo com 5 pontos (até $10^6$), obtemos **$c_0 = -0.5298$** e **$c_1 = 0.9533$**, que são então avaliados no ponto de teste $x = 10^7$ (Seção 3).

---

### 2.4 Aspectos Computacionais: Chave RSA de 1024 bits ($x = 2^{1024}$) e Gaveta de 64 bits vs. `BigFloat`

Uma chave RSA de 1024 bits possui 1024 posições binárias (`0` ou `1`), abrangendo inteiros até $x = 2^{1024}$.
- **Estouro (Overflow) em 64 bits:** Tanto o tipo `number` do TypeScript quanto o `Float64` da Julia utilizam registros fixos de **64 bits** (1 bit de sinal, 11 bits de expoente, 52 bits de mantissa), cujo maior valor representável é $\approx 2^{1023} \approx 1.797 \times 10^{308}$. Ao executar `2 ** 1024` diretamente em TypeScript, ocorre **Overflow** e a linguagem retorna `Infinity`.
- **Resolução Direta na Julia (`BigFloat`):** A Julia dispõe nativamente do tipo `BigFloat` (alocado dinamicamente em RAM via biblioteca GNU MPFR). Executando `BigFloat(2)^1024`, a Julia calcula o valor sem overflow ($1.7977 \times 10^{308}$) e obtém diretamente $\ln(2^{1024}) = 709.7827$.
- **Resolução Algébrica no TypeScript:** Para contornar o `Infinity` sem biblioteca externa, aplicamos a propriedade de logaritmos antes de avaliar a potência:

$$\ln(2^{1024}) = 1024 \cdot \ln(2) = 1024 \times 0.693147 = 709.7827$$

Com $\ln(x) = 709.7827$, respondemos à pergunta-gancho do projeto:
1. **Densidade de primos perto de $2^{1024}$:** Pela relação $\frac{1}{\ln(x)} \approx \frac{1}{709.78}$, cerca de **1 em cada 710 números** de 1024 bits é primo ($0.141\%$).
2. **Tentativas médias do RSA:** Como o gerador de chaves descarta números pares de imediato, ele precisa testar em média **$709.78 / 2 \approx 355$ candidatos ímpares** até encontrar um primo.
3. **Total de primos de 1024 bits (entre $2^{1023}$ e $2^{1024}$):** Existem aproximadamente $\frac{2^{1023}}{709.78} \approx \mathbf{1.266 \times 10^{305}}$ números primos de 1024 bits, garantindo que diferentes chaves RSA geradas aleatoriamente jamais colidam.

---

## 3. Discussão e Conclusão

### 3.1 Comparação dos Modelos na Escala Original de $\pi(x)$ (Determinístico)

Um cuidado fundamental na linearização (Aula 04) é que o Sistema Normal transformado minimiza os desvios na variável $w = x/\pi(x)$, e não diretamente em $\pi(x)$. Portanto, para comparar modelos de maneira justa, todas as previsões devem ser convertidas de volta para a **escala original de $\pi(x)$**. A Tabela 4 resume o desempenho de cada aproximação no ponto de teste $x = 10^7$ (que ficou fora do conjunto de treino $10^2..10^6$).

**Tabela 4 — Comparação das aproximações no ponto de teste $x = 10^7$ ($\pi(10^7)$ real = $664\,579$):**

| Modelo / Aproximação | Fórmula / Parâmetros | $\pi(10^7)$ Estimado | Erro Absoluto | Erro Relativo (%) |
|---|---|---|---|---|
| Assintótico puro (TNP) | $x / \ln(x)$ | $620\,421$ | $-44\,158$ | $-6.64\%$ |
| Legendre Linearizado (Treino $10^2..10^6$) | $c_0 = -0.5298,\; c_1 = 0.9533$ | $674\,031$ | $+9\,452$ | $+1.42\%$ |
| Mínimos Quadrados Linear (Treino $10^2..10^6$) | $a = -53.40,\; b = 0.4706$ | $672\,531$ | $+7\,952$ | $+1.20\%$ |
| Legendre Histórico (1808) | $c_0 = -1.08366,\; c_1 = 1.0000$ | $665\,140$ | $+561$ | **$+0.08\%$** |
| **Gabarito Real (Crivo)** | **Contagem exata** | **$664\,579$** | — | — |

**Por que os nossos parâmetros ajustados ($c_0 = -0.5298, c_1 = 0.9533$) erraram $+1.42\%$ em $10^7$, enquanto as constantes originais de Legendre ($c_0 = -1.08366, c_1 = 1$) erraram apenas $+0.08\%$?**
Por duas razões matemáticas:
1. O nosso ajuste incluiu pontos pequenos ($x = 10^2$ e $10^3$), onde a função $\pi(x)$ ainda não entrou no regime assintótico suave, "puxando" a inclinação $c_1$ para baixo de $1$.
2. A linearização atribui o mesmo peso quadrático aos desvios em $x/\pi(x)$ tanto em $x = 100$ quanto em $x = 10^6$, em vez de minimizar o erro relativo em $\pi(x)$ para $x$ grande.

---

### 3.2 Benchmark Pareado Não-Determinístico: TypeScript (V8) vs. Julia (LLVM)

Como o tempo de execução sofre variações não-determinísticas do sistema operacional e do coletor de lixo (gerando distribuições assimétricas com cauda à direita), aplicamos um protocolo pareado idêntico nas duas linguagens:
- **Aquecimento obrigatório:** **5 execuções iniciais descartadas** (eliminando a latência de compilação JIT do LLVM na Julia e do TurboFan no V8 do Node.js).
- **Amostragem:** **30 execuções consecutivas medidas**, reportando a **Mediana**, o **Mínimo** e a **Média**, além do consumo de memória da estrutura de dados.

**Tabela 5 — Resultados do Benchmark Pareado (5 warmups descartados + 30 amostras, Apple Silicon):**

| Operação Avaliada | Métrica | TypeScript (Node v24 / V8) | Julia (v1.13.3 / LLVM) | Análise Técnica |
|---|---|---|---|---|
| **Crivo de Eratóstenes ($10^7$)** | Estrutura e Memória | `Uint8Array` (**9.54 MB**) | `BitVector` (**1.19 MB**) | Julia usa **8× menos memória** (1 bit vs. 1 byte por número) |
| **Crivo de Eratóstenes ($10^7$)** | Tempo Mediana (30 runs) | **38.70 ms** | **43.26 ms** | Trade-off Tempo × Espaço (ver nota abaixo) |
| **Crivo de Eratóstenes ($10^7$)** | Tempo Mínimo (30 runs) | **36.54 ms** | **42.60 ms** | Ambas na faixa de ~40 ms após aquecimento |
| **Mínimos Quadrados (Ajustes)** | Método de Resolução | Sistema Normal + Gauss manual | Operador nativo `S \ y` | Julia resolve em 1 linha sem implementar Gauss |
| **Mínimos Quadrados (Ajustes)** | Tempo Mediana (30 runs) | **9.85 µs** | **3.92 µs** | Julia é **2.5× mais rápida** na mediana |
| **Mínimos Quadrados (Ajustes)** | Tempo Mínimo (30 runs) | **9.37 µs** | **1.71 µs** | Julia é **5.5× mais rápida** no tempo mínimo |

**Discussão do Trade-off Tempo × Espaço no Crivo:**
Os dados da Tabela 5 revelam um clássico trade-off de Ciência da Computação: no Crivo até $10^7$, o `Uint8Array` do TypeScript teve mediana ligeiramente menor ($38.70\text{ ms}$ contra $43.26\text{ ms}$ da Julia), porém consumiu **8 vezes mais memória RAM ($9.54\text{ MB}$ contra $1.19\text{ MB}$)**. Isso ocorre porque o `Uint8Array` grava 1 byte inteiro direto na memória por instrução, enquanto o `BitVector` da Julia compacta 64 números em cada palavra de 64 bits, exigindo operações extras de deslocamento de bits (*bit-shifting* e máscara booleana) a cada marcação em troca de economizar 87.5% de memória. Já na resolução matricial dos Mínimos Quadrados, as rotinas nativas de Álgebra Linear da Julia (`S \ y`) superaram a implementação manual em TypeScript em **2.5× na mediana ($3.92\text{ µs}$ vs. $9.85\text{ µs}$)** e **5.5× no tempo mínimo ($1.71\text{ µs}$ vs. $9.37\text{ µs}$)**.

### 3.3 Conclusão

O Método dos Mínimos Quadrados, combinado à técnica de linearização, permitiu reconstruir a partir de apenas 4 a 5 pontos empíricos do Crivo de Eratóstenes aproximações para $\pi(x)$ com erro inferior a $2\%$ em $10^6$ e $10^7$. O comparativo entre TypeScript e Julia evidenciou que, enquanto linguagens de propósito geral exigem atenção manual do desenvolvedor contra armadilhas do padrão IEEE 754 de 64 bits (como cancelamento subtrativo em matrizes mal-condicionadas e *overflow* em $2^{1024}$), linguagens científicas oferecem abstrações nativas (`S \ y` e `BigFloat`) que garantem concisão, velocidade matricial e precisão arbitrária.

---

## 4. Bibliografia

- BEZANSON, J.; EDELMAN, A.; KARPINSKI, S.; SHAH, V. B. Julia: A fresh approach to numerical computing. *SIAM Review*, v. 59, n. 1, p. 65–98, 2017.
- BURDEN, R. L.; FAIRES, J. D.; BURDEN, A. M. *Análise Numérica*. 10. ed. São Paulo: Cengage Learning, 2016.
- CORMEN, T. H.; LEISERSON, C. E.; RIVEST, R. L.; STEIN, C. *Algoritmos: Teoria e Prática*. 3. ed. Rio de Janeiro: Elsevier, 2012. (Cap. 31: Algoritmos de Teoria dos Números e o Sistema RSA.)
- FRANCO, N. B. *Cálculo Numérico*. São Paulo: Pearson Prentice Hall, 2006. (Cap. 1: Erros; Cap. 3: Sistemas Lineares; Cap. 7: Aproximação de Funções — Método dos Mínimos Quadrados, Problema 7.3.)
- MORAIS, H. *Notas de aula de Cálculo Numérico (Aulas 01, 03 e 04)*. Rio Claro: IGCE/DEMAC — UNESP, 2026.
