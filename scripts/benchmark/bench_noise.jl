include("bench_utils.jl");

# =============================================================================
# Benchmark BCFEigen construction.
# =============================================================================

let
    title = "Benchmark BCFEigen construction."
    print_header(title)

    # Construct BCF.
    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)

    @info "Three-mode squeezed reservoir BCF."

    t_end     = 20.0
    grid_size = 1000

    bcf_eigen = Humulus.BCFEigen(bcf, t_end, grid_size)

    @info "Benchmark parameters:" grid_size

    bench = @benchmark Humulus.BCFEigen($bcf, $t_end, $grid_size)
    benchmark_construction(bcf_eigen, bench)
end


# =============================================================================
# Benchmark NoiseSampler construction.
# =============================================================================

let
    title = "Benchmark NoiseSampler construction."
    print_header(title)

    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)

    @info "Three-mode squeezed reservoir BCF."

    t_end = 20.0
    grid_size = 5_000

    bcf_eigen = Humulus.BCFEigen(bcf, t_end, grid_size)

    mkpath("bcf_eigen_cache")
    path = joinpath("bcf_eigen_cache", "bench_noise.jld2")
    save_object(path, bcf_eigen)

    sampler! = Humulus.sampler_from_cache(path; checks = true, clear_cache = false);

    @info "Benchmark parameters" grid_size

    bench = @benchmark Humulus.sampler_from_cache($path; checks = true, clear_cache = false)
    benchmark_construction(sampler!, bench)
    rm(path)
end


# =============================================================================
# Benchmark NoiseSampler evaluation.
# =============================================================================

let
    title = "Benchmark NoiseSampler evaluation."
    print_header(title)

    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)

    @info "Three-mode squeezed reservoir BCF."

    sampler! = let
        t_end = 20.0
        grid_size = 500

        bcf_eigen = Humulus.BCFEigen(bcf, t_end, grid_size)

        mkpath("bcf_eigen_cache")
        path = joinpath("bcf_eigen_cache", "bench_noise.jld2")

        save_object(path, bcf_eigen)

        Humulus.sampler_from_cache(path, checks = true, clear_cache = true)
    end

    grid_size = sampler!._time_grid.n_points

    @info "Benchmark parameters" grid_size

    noise = zeros(ComplexF64, grid_size)

    bench = @benchmark $sampler!($noise)
    benchmark_evaluation(bench)
end