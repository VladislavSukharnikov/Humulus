testset_name = "HOPS";
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

        noise_params = NoiseParams(t_end, bcf, grid_size);

        mkpath("noise_cache")
        path = joinpath("noise_cache", "test_hops.jld2")
        save_object(path, noise_params)

        # Construct SampleNoise object.
        sample_noise! = load_noise(path);
        rm(path)
        sample_noise!
    end

    # Number of effective modes.
    N = n_modes(bcf)

    # Maximum occupation number allowed for each mode.
    max_occupancies = (5,5,5)

    # HOPS construction.
    max_fock_states = maximum(max_occupancies)+1
    hops! = HOPS{N,max_fock_states}(sample_noise!._time_grid)
    
    @test hops! isa HOPS{N, max_fock_states}

    println("   $testset_name: testing type parameters...")
    @testset "Type parameters" begin
        @test typeof(hops!).parameters[1] == N
        @test typeof(hops!).parameters[2] == max_fock_states
    end

    println("   $testset_name: testing fields...")
    @testset "Fields" begin
        for idx in 1:max_fock_states
            @test hops!._sqrt_of_int[idx]==sqrt(idx-1)
        end
        @test hops!.f_tmp isa MVector{N, ComplexF64}
        @test hops!.g_tmp isa MVector{N, ComplexF64}
        @test hops!._time_grid == sample_noise!._time_grid
        @test hops!.noise isa Vector{ComplexF64}
        @test length(hops!.noise) == length(hops!._time_grid)
    end

    println("  $testset_name: test complete.")
end;

testset_name = "HOPS (function)"
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

        noise_params = NoiseParams(t_end, bcf, grid_size);

        mkpath("noise_cache")
        path = joinpath("noise_cache", "test_hops.jld2")
        save_object(path, noise_params)

        # Construct SampleNoise object.
        sample_noise! = load_noise(path);
        rm(path)
        sample_noise!
    end

    atom_params   = let
        ν₀            = rand()
        c             = rand(ComplexF64, 2)
        c           ./= sqrt(abs2(c[1])+abs2(c[2]))
        atom_params   = AtomParams(ν₀, c...)
    end

    # Number of effective modes.
    N = n_modes(bcf)

    # Maximum occupation number allowed for each mode.
    max_occupancies = (5,5,5)

    # FockSpace construction.
    fock_space = FockSpace(Val(N), max_occupancies)

    # HOPS construction.
    max_fock_states = maximum(max_occupancies)+1
    hops! = HOPS{N,max_fock_states}(sample_noise!._time_grid)
    
    # Check that all elements of noise before sampling is zero.
    @test all(iszero, hops!.noise)

    # Sample noise. 
    sample_noise!(hops!.noise)

    # Check that all elements of noise are not zero now.
    @test all(!iszero, hops!.noise)


    fock_dim = fock_space.fock_dim

    dψ = rand(ComplexF64, 2, fock_dim)
    ψ  = rand(ComplexF64, 2, fock_dim)

    dX = rand(ComplexF64, N, 1)
    X  = rand(ComplexF64, N, 1)

    du = ArrayPartition(dψ, dX)
    u  = ArrayPartition(ψ, X)

    solver_params = create_solver_params(bcf, fock_space, atom_params)
    t = rand()*sample_noise!._time_grid[end]


    println("   $testset_name: testing type stability...")
    @testset "Type stability" begin
        @inferred hops!(du, u, solver_params, t)
    end


    println("   $testset_name: testing allocations...")
    @testset "Allocations" begin
        hops!(du, u, solver_params, t)
        bench = @benchmark $hops!($du, $u, $solver_params, $t)
        @test mean(bench).memory==0
    end

    println("   $testset_name: testing mutation...")
    @testset "Mutation test" begin
        @. du  = 0.0

        u_copy  = deepcopy(u)
        du_copy = deepcopy(du)
        solver_params_copy = deepcopy(solver_params)

        hops!(du_copy, u_copy, solver_params_copy, t)

        @test u_copy==u
        @test solver_params==solver_params_copy
        @test du_copy != du

        hops!(du, u_copy, solver_params_copy, t)

        @test du == du_copy
    end
    
    println("  $testset_name: test complete.")
end;
