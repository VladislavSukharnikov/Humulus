@testset "BCFEigen" begin
    @info "BCFEigen: construction..."

    test_cases = (
        (1, 100),
        (2, 150),
        (3, 200),
    )

    @testset "N=$N, grid_size=$grid_size" for (N, grid_size) in test_cases
        bcf = random_bcf(N)

        t_end = 20.0
        bcf_eigen = BCFEigen(bcf, t_end, grid_size)

        mkpath("covariance_cache")
        path = joinpath("covariance_cache", "test.jld2")
        save_object(path, bcf_eigen)

        @testset "Serialization" begin
            loaded_bcf_eigen = load_object(path)

            @test loaded_bcf_eigen.bcf == bcf_eigen.bcf
            @test loaded_bcf_eigen.time_grid == bcf_eigen.time_grid
            @test loaded_bcf_eigen.vecs == bcf_eigen.vecs
        end

        rm(path; force=true)
    end
    @info "BCFEigen: completed."
end;


@testset "BCFCholesky" begin
    @info "BCFCholesky: construction..."

    test_cases = (
        (1, 100),
        (2, 150),
        (3, 200),
    )

    @testset "N=$N, grid_size=$grid_size" for (N, grid_size) in test_cases
        bcf = random_bcf(N)

        t_end = 20.0
        bcf_chol = BCFCholesky(bcf, t_end, grid_size)

        mkpath("covariance_cache")
        path = joinpath("covariance_cache", "test.jld2")
        save_object(path, bcf_chol)

        @testset "Serialization" begin
            loaded_bcf_chol = load_object(path)

            @test loaded_bcf_chol.bcf == bcf_chol.bcf
            @test loaded_bcf_chol.time_grid == bcf_chol.time_grid
            @test loaded_bcf_chol.chol == bcf_chol.chol
        end

        rm(path; force=true)
    end
    @info "BCFCholesky: completed."
end;


@testset "NoiseSampler" begin
    @info "NoiseSampler: construction..."

    test_cases = (
        (1, 100),
        (2, 150),
        (3, 200),
    )

    for Decomp in (BCFEigen, BCFCholesky)

        @testset "$(Decomp), N=$N, grid_size=$grid_size" for (N, grid_size) in test_cases

            bcf = random_bcf(N)

            t_end = 20.0
            decomposition = Decomp(bcf, t_end, grid_size)

            mkpath("covariance_cache")
            path = joinpath("covariance_cache", "test.jld2")
            save_object(path, decomposition)

            sampler! = sampler_from_cache(path; checks=true, clear_cache=true)

            @test sampler!._time_grid == decomposition.time_grid
            @test sampler!._A == _noise_factor(decomposition, true)

        end
    end
end;


@testset "NoiseSampler (function)" begin
    @info "NoiseSampler: evaluation..."

    test_cases = (
        (1, 100),
        (2, 150),
        (3, 200),
    )

    for Decomp in (BCFEigen, BCFCholesky)

        @testset "$(Decomp), N=$N, grid_size=$grid_size" for (N, grid_size) in test_cases

            bcf = random_bcf(N)

            t_end = 20.0
            decomposition = Decomp(bcf, t_end, grid_size)

            mkpath("covariance_cache")
            path = joinpath("covariance_cache", "test.jld2")
            save_object(path, decomposition)

            sampler! = sampler_from_cache(path; checks=true, clear_cache=true)

            @testset "Type stability" begin
                noise = zeros(ComplexF64, grid_size)
                @inferred sampler!(noise)
            end

            @testset "Sampling" begin
                noise = zeros(ComplexF64, grid_size)

                sampler!(noise)

                @test any(!iszero, noise)
                @test all(isfinite, real.(noise))
                @test all(isfinite, imag.(noise))
                @test all(!isnan, real.(noise))
                @test all(!isnan, imag.(noise))
            end
        end

        @testset "Allocations ($(Decomp))" begin
            N = 1
            grid_size = 100

            bcf = random_bcf(N)

            t_end = 20.0
            decomposition = Decomp(bcf, t_end, grid_size)

            mkpath("covariance_cache")
            path = joinpath("covariance_cache", "test.jld2")
            save_object(path, decomposition)

            sampler! = sampler_from_cache(path; checks=true, clear_cache=true)

            noise = zeros(ComplexF64, grid_size)

            sampler!(noise)

            @test @ballocated($sampler!($noise)) == 0
        end

    end

    @info "NoiseSampler: completed."
end;