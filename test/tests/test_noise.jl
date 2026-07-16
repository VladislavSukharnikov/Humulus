@testset "BCFEigen" begin
    @info "BCFEigen: construction"

    test_cases = (
        (1, 100),
        (2, 150),
        (3, 200),
    )

    @testset "N=$N, grid_size=$grid_size" for (N, grid_size) in test_cases
        bcf = random_bcf(N)

        t_end = 20.0
        bcf_eigen = BCFEigen(bcf, t_end, grid_size)

        mkpath("bcf_eigen_cache")
        path = joinpath("bcf_eigen_cache", "test.jld2")
        save_object(path, bcf_eigen)

        @testset "Serialization" begin
            loaded_bcf_eigen = load_object(path)

            @test loaded_bcf_eigen.bcf == bcf_eigen.bcf
            @test loaded_bcf_eigen.time_grid == bcf_eigen.time_grid
            @test loaded_bcf_eigen.vals == bcf_eigen.vals
            @test loaded_bcf_eigen.vecs == bcf_eigen.vecs
        end

        rm(path; force=true)
    end
    @info "BCFEigen: done"
end;


@testset "NoiseSampler" begin
    @info "NoiseSampler: construction"

    test_cases = (
        (1, 100),
        (2, 150),
        (3, 200),
    )

    @testset "N=$N, grid_size=$grid_size" for (N, grid_size) in test_cases
        bcf = random_bcf(N)

        t_end = 20.0
        bcf_eigen = BCFEigen(bcf, t_end, grid_size)

        mkpath("bcf_eigen_cache")
        path = joinpath("bcf_eigen_cache", "test.jld2")
        save_object(path, bcf_eigen)

        @testset "Construction" begin
            sampler! = sampler_from_cache(path; checks=true, clear_cache=true)

            @test sampler!._time_grid == bcf_eigen.time_grid
            @test sampler!._vals == bcf_eigen.vals
            @test sampler!._vecs == bcf_eigen.vecs
        end
    end
end;


@testset "NoiseSampler (function)" begin
    @info "NoiseSampler: evaluation"

    test_cases = (
        (1, 100),
        (2, 150),
        (3, 200),
    )

    @testset "N=$N, grid_size=$grid_size" for (N, grid_size) in test_cases
        bcf = random_bcf(N)

        t_end = 20.0
        bcf_eigen = BCFEigen(bcf, t_end, grid_size)

        mkpath("bcf_eigen_cache")
        path = joinpath("bcf_eigen_cache", "test.jld2")
        save_object(path, bcf_eigen)

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

    @testset "Allocations" begin
        N = 1
        grid_size = 100

        bcf = random_bcf(N)

        t_end = 20.0
        bcf_eigen = BCFEigen(bcf, t_end, grid_size)

        mkpath("bcf_eigen_cache")
        path = joinpath("bcf_eigen_cache", "test.jld2")
        save_object(path, bcf_eigen)

        sampler! = sampler_from_cache(path; checks=true, clear_cache=true)

        noise = zeros(ComplexF64, grid_size)

        sampler!(noise)

        @test @ballocated($sampler!($noise)) == 0
    end
    @info "NoiseSampler: done"
end;