include("bench_utils.jl");

# =============================================================================
# Benchmark random BCF construction.
# =============================================================================

let
    N = 100

    title = "Construction memory of the random $N-mode BCF."
    print_header(title)


    @info "Benchmark parameters:" N

    bcf   = Humulus.random_bcf(N)

    bench = @benchmark Humulus.random_bcf($N)
    benchmark_construction(bcf, bench)
end

# =============================================================================
# Benchmark random BCF evaluation.
# =============================================================================

let
    N = 100

    title = "Benchmark random $N-mode BCF evaluation."
    print_header(title)


    @info "Benchmark parameters:" N

    bcf = Humulus.random_bcf(N)
    t, s = rand(), rand()
    bcf(t, s)

    bench = @benchmark $bcf($t, $s)
    benchmark_evaluation(bench)
end


# =============================================================================
# Benchmark one-mode squeezed-reservoir BCF evaluation.
# =============================================================================

let
    title = "Benchmark single-mode squeezed-reservoir BCF evaluation."
    print_header(title)

    bcf = one_mode_squeezed_bcf(5.0, 1.0, 1.5, 0.0, 1.0)
    t,s = rand(),rand()
    bcf(t, s)

    bench = @benchmark $bcf($t, $s)
    benchmark_evaluation(bench)
end


# =============================================================================
# Benchmark three-mode squeezed-reservoir BCF evaluation.
# =============================================================================

let
    title = "Benchmark three-mode squeezed-reservoir BCF evaluation."
    print_header(title)

    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)
    t,s = rand(),rand()
    bcf(t, s)

    bench = @benchmark $bcf($t, $s)
    benchmark_evaluation(bench)
end