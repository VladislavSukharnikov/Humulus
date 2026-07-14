testset_name = "NoiseParams";
@testset "$testset_name" begin
    println("   $testset_name: setting сonfiguration...")

    bcf = let 
        ω₀ = 5.0;         # central (carrier) frequency
        Γ₀ = 2.0;         # input spectral width
        κ  = 1.0;         # cavity loss rate
        ϵ  = 0.5;         # effective DPA pump amplitude
        φ  = Float64(π);  # squeezing phase
        γ  = 1.0;         # atom–field coupling strength

        three_mode_squeezed_bcf(ω₀, Γ₀, κ, ϵ, φ, γ);
    end

    t_end        = 20.0
    grid_size    = 200
    noise_params = NoiseParams(bcf, t_end, grid_size)

    mkpath("noise_cache")
    path = joinpath("noise_cache", "test_$testset_name.jld2")
    save_object(path, noise_params)

    println("   $testset_name: testing loading...")
    @testset "Load NoiseParams" begin
        loaded_noise_params = load_object(path)
        @test loaded_noise_params.bcf == noise_params.bcf
        @test loaded_noise_params.vals == noise_params.vals
        @test loaded_noise_params.vecs == noise_params.vecs
        @test loaded_noise_params.time_grid == noise_params.time_grid
        @test loaded_noise_params.grid_size == noise_params.grid_size
    end

    rm(path; force = true)

    println("   $testset_name: test complete.")
end;

testset_name = "SampleNoise";
@testset "$testset_name" begin
    println("   $testset_name: setting сonfiguration...")
    
    bcf = let 
        ω₀ = 5.0;         # central (carrier) frequency
        Γ₀ = 2.0;         # input spectral width
        κ  = 1.0;         # cavity loss rate
        ϵ  = 0.5;         # effective DPA pump amplitude
        φ  = Float64(π);  # squeezing phase
        γ  = 1.0;         # atom–field coupling strength

        three_mode_squeezed_bcf(ω₀, Γ₀, κ, ϵ, φ, γ);
    end

    t_end        = 20.0
    grid_size    = 200
    noise_params = NoiseParams(bcf, t_end, grid_size)

    mkpath("noise_cache")
    path = joinpath("noise_cache", "test_$testset_name.jld2")
    save_object(path, noise_params)

    println("   $testset_name: testing construction...")
    @testset "SampleNoise construction" begin
        sample_noise! = load_noise(path)
        @test sample_noise!._time_grid == noise_params.time_grid
        @test sample_noise!._grid_size == noise_params.grid_size
        @test sample_noise!._vals == noise_params.vals
        @test sample_noise!._vecs == noise_params.vecs
        @test sample_noise!._path == path
    end

    rm(path; force = true)

    println("   $testset_name: test complete.")
end

testset_name = "SampleNoise (function)";
@testset "$testset_name" begin
    println("   $testset_name: setting сonfiguration...")
    
    bcf = let 
        ω₀ = 5.0;         # central (carrier) frequency
        Γ₀ = 2.0;         # input spectral width
        κ  = 1.0;         # cavity loss rate
        ϵ  = 0.5;         # effective DPA pump amplitude
        φ  = Float64(π);  # squeezing phase
        γ  = 1.0;         # atom–field coupling strength

        three_mode_squeezed_bcf(ω₀, Γ₀, κ, ϵ, φ, γ);
    end

    sample_noise! = let 
        
        t_end     = 20.0;
        grid_size = 200;

        noise_params = NoiseParams(bcf, t_end, grid_size);

        mkpath("noise_cache")
        path = joinpath("noise_cache", "test_$testset_name.jld2")
        save_object(path, noise_params)

        # Construct SampleNoise object.
        sample_noise! = load_noise(path);
        rm(path)
        sample_noise!
    end

    println("   $testset_name: testing type stability...")
    @testset "Type stability" begin
        z = zeros(ComplexF64, sample_noise!._grid_size)
        @inferred sample_noise!(z)
    end

    println("   $testset_name: testing allocations...")
    @testset "Allocations" begin
        z = zeros(ComplexF64, sample_noise!._grid_size)

        sample_noise!(z)
        bench = @benchmark $sample_noise!($z)
        @test mean(bench.memory) == 0
    end

    println("   $testset_name: test complete.")
end