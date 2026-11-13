using LinearAlgebra
using Printf
using Statistics

println("====================================================================")
println(" PROJETO CÁLCULO NUMÉRICO — ESTIMANDO π(x) E CHAVES RSA (JULIA)     ")
println("====================================================================\n")

function sieve_primes(limit::Int, checkpoints::Vector{Int})
    is_prime = trues(limit)
    is_prime[1] = false

    sqrt_limit = isqrt(limit)
    @inbounds for p in 2:sqrt_limit
        if is_prime[p]
            for mult in (p * p):p:limit
                is_prime[mult] = false
            end
        end
    end

    counts = Dict{Int, Int}()
    running = 0
    cp_idx = 1
    sorted_cps = sort(checkpoints)

    @inbounds for i in 1:limit
        if is_prime[i]
            running += 1
        end
        while cp_idx <= length(sorted_cps) && i == sorted_cps[cp_idx]
            counts[sorted_cps[cp_idx]] = running
            cp_idx += 1
        end
    end

    bytes_used = sizeof(is_prime.chunks)
    return counts, bytes_used
end

function run_paired_benchmark(fn::Function, warmups::Int=5, runs::Int=30)
    for _ in 1:warmups
        fn()
    end
    samples_ms = Vector{Float64}(undef, runs)
    for i in 1:runs
        t0 = time_ns()
        fn()
        samples_ms[i] = (time_ns() - t0) / 1e6
    end
    sort!(samples_ms)
    return (
        min_ms = samples_ms[1],
        median_ms = median(samples_ms),
        mean_ms = mean(samples_ms)
    )
end

checkpoints = [10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000]
counts, bytes_used = sieve_primes(10_000_000, checkpoints)

# ---------------------------------------------------------------------------
# BLOCO 1: MÍNIMOS QUADRADOS LINEAR E LINEARIZAÇÃO DE LEGENDRE (S \ y)
# ---------------------------------------------------------------------------
x_book = Float64[100, 1_000, 10_000, 100_000]
y_book = Float64[counts[Int(x)] for x in x_book]
u_book = [x / log10(x) for x in x_book]
S_book = hcat(ones(length(u_book)), u_book)
a_lin, b_lin = S_book \ y_book
residuo_Q = sum(((S_book * [a_lin, b_lin]) .- y_book) .^ 2)
pred_6dig = (a_lin + b_lin * (1e6 / 6.0)) - (a_lin + b_lin * (1e5 / 5.0))
real_6dig = counts[1_000_000] - counts[100_000]

println("[1. PONTE 1] Mínimos Quadrados Linear (Operador nativo S \\ y)")
@printf("  • Equação ajustada  : π(x) ≈ %.2f + %.4f * (x / log10(x)) | Q = %.1f\n", a_lin, b_lin, residuo_Q)
@printf("  • Primos 6 dígitos  : %d estimados vs. %d reais (+1.97%%)\n\n", round(Int, pred_6dig), real_6dig)

x_train = Float64[100, 1_000, 10_000, 100_000, 1_000_000]
y_train = Float64[counts[Int(x)] for x in x_train]
c0_4, c1_4 = hcat(ones(length(x_book)), log.(x_book)) \ (x_book ./ y_book)
c0_tr, c1_tr = hcat(ones(length(x_train)), log.(x_train)) \ (x_train ./ y_train)
pred_1e6_4pts = 1e6 / (c0_4 + c1_4 * log(1e6))
pred_1e7_tr = 1e7 / (c0_tr + c1_tr * log(1e7))
pred_1e7_1808 = 1e7 / (log(1e7) - 1.08366)

println("[2. PONTE 2] Linearização de Legendre: w = x/π(x) = c0 + c1*ln(x)")
@printf("  • 4 pontos (10²..10⁵): c0 = %.4f, c1 = %.4f -> π(10⁶) ≈ %d (real: 78498, +1.96%%)\n",
        c0_4, c1_4, round(Int, pred_1e6_4pts))
@printf("  • Treino até 10⁶     : c0 = %.4f, c1 = %.4f -> π(10⁷) ≈ %d (real: 664579, +1.42%%)\n",
        c0_tr, c1_tr, round(Int, pred_1e7_tr))
@printf("  • Legendre (1808)    : c0 = -1.08366, c1 = 1.0000 -> π(10⁷) ≈ %d (real: 664579, +0.08%%)\n\n",
        round(Int, pred_1e7_1808))

# ---------------------------------------------------------------------------
# BLOCO 2 (PARA PRINT NO SLIDE): RSA 1024 BITS — Float64 vs. BigFloat
# ---------------------------------------------------------------------------
float64_2_1024 = Float64(2)^1024
setprecision(BigFloat, 256)
big_2_1024 = BigFloat(2)^1024
ln_big = log(big_2_1024)
primes_1024_bits = (BigFloat(2)^1023) / ln_big

println("[3. RSA 1024 BITS] Gaveta Fixa Float64 vs. Gaveta Elástica BigFloat")
@printf("  • Float64(2)^1024 (64 bits IEEE 754) : %s (Overflow!)\n", string(float64_2_1024))
@printf("  • BigFloat(2)^1024 (256 bits na RAM) : %.4e (Calculado direto sem overflow!)\n", big_2_1024)
@printf("  • Sorteios ímpares médios no RSA     : %d tentativas (1 primo a cada %d números)\n",
        round(Int, ln_big / 2), round(Int, ln_big))
@printf("  • Total de primos de 1024 bits       : %.3e primos\n\n", primes_1024_bits)

# ---------------------------------------------------------------------------
# BLOCO 3 (PARA PRINT NO SLIDE): BENCHMARK PAREADO (5 WARMUPS + 30 EXECUÇÕES)
# ---------------------------------------------------------------------------
bench_sieve = run_paired_benchmark(() -> sieve_primes(10_000_000, checkpoints), 5, 30)

function fit_both_julia(x_tr::Vector{Float64}, y_tr::Vector{Float64})
    u_tr = x_tr ./ log10.(x_tr)
    p_lin = hcat(ones(length(u_tr)), u_tr) \ y_tr
    X_tr = log.(x_tr)
    w_tr = x_tr ./ y_tr
    p_leg = hcat(ones(length(X_tr)), X_tr) \ w_tr
    return p_lin, p_leg
end

bench_ls = run_paired_benchmark(() -> fit_both_julia(x_train, y_train), 5, 30)

println("[4. BENCHMARK JULIA] Performance Não-Determinística (5 warmups + 30 amostras)")
@printf("  • Crivo 10^7 (BitVector %.2f MB) : Mediana = %.2f ms | Mínimo = %.2f ms | Média = %.2f ms\n",
        bytes_used / (1024 * 1024), bench_sieve.median_ms, bench_sieve.min_ms, bench_sieve.mean_ms)
@printf("  • Mínimos Quadrados (S \\ y)     : Mediana = %.2f µs | Mínimo = %.2f µs\n",
        bench_ls.median_ms * 1000, bench_ls.min_ms * 1000)
