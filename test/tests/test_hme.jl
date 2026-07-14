testset_name = "HME"
@testset "$testset_name" begin
    println("   $testset_name: setting configuration...")

    # Number of effective modes.
    N = 2

    max_occupancies = (10, 10)
    max_fock_states = maximum(max_occupancies) + 1

    fock_space = FockSpace(Val(N), max_occupancies)
    bcf           = random_bcf(N)
    atom_params   = let
        ν₀            = rand()
        c             = rand(ComplexF64, 2)
        c           ./= sqrt(abs2(c[1])+abs2(c[2]))
        atom_params   = AtomParams(ν₀, c...)
    end
    solver_params = create_solver_params(bcf, fock_space, atom_params)
    
    hme! = HME{N, max_fock_states}()
    @test hme! isa HME{N, max_fock_states}

    println("   $testset_name: testing type parameters...")
    @testset "Type parameters" begin
        @test typeof(hme!).parameters[1] == N
        @test typeof(hme!).parameters[2] == max_fock_states
    end

    println("   $testset_name: testing fields...")
    @testset "Fields" begin
        for idx in 1:max_fock_states
            @test hme!._sqrt_of_int[idx]==sqrt(idx-1)
        end
        @test hme!.f_tmp isa MVector{N, ComplexF64}
        @test hme!.g_tmp isa MVector{N, ComplexF64}
    end

    println("   $testset_name: test complete.")
end;

testset_name = "HME (function)";
@testset "$testset_name" begin
    println("   $testset_name: setting configuration...")

    # Number of effective modes.
    N = 2

    max_occupancies = (10, 10)
    max_fock_states = maximum(max_occupancies) + 1

    fock_space = FockSpace(Val(N), max_occupancies)
    bcf           = random_bcf(N)
    atom_params   = let
        ν₀            = rand()
        c             = rand(ComplexF64, 2)
        c           ./= sqrt(abs2(c[1])+abs2(c[2]))
        atom_params   = AtomParams(ν₀, c...)
    end
    solver_params = create_solver_params(bcf, fock_space, atom_params)
    
    hme! = HME{N, max_fock_states}()

    fock_dim = solver_params.fock_dim
    t        = rand()
    ρ        = rand(ComplexF64, 2, 2, fock_dim, fock_dim)
    dρ       = rand(ComplexF64, 2, 2, fock_dim, fock_dim)

    println("   $testset_name: testing type stability...")
    @testset "Type stability" begin
        @inferred hme!(dρ, ρ, solver_params, t)
    end

    println("   $testset_name: testing evaluation...")
    @testset "Evaluation" begin
        hme!(dρ, ρ, solver_params, t)
        @test all(isfinite, real.(dρ))
        @test all(isfinite, imag.(dρ))
    end

    println("   $testset_name: testing allocations...")
    @testset "Allocations" begin
        hme!(dρ, ρ, solver_params, t)

        bench = @benchmark $hme!($dρ, $ρ, $solver_params, $t)
        @test mean(bench).memory==0
    end

    println("   $testset_name: testing mutation...")
    @testset "Mutation test" begin
        @. dρ  = 0.0

        ρ_copy  = deepcopy(ρ)
        dρ_copy = deepcopy(dρ)
        solver_params_copy = deepcopy(solver_params)

        hme!(dρ_copy, ρ_copy, solver_params_copy, t)

        @test ρ_copy==ρ
        @test solver_params==solver_params_copy
        @test dρ_copy != dρ

        hme!(dρ, ρ_copy, solver_params_copy, t)

        @test dρ == dρ_copy
    end

    println("   $testset_name: test complete.")
end;