@testset "HME" begin
    println("HME: construction")

    test_cases = (
        (1, 0),
        (1, (5,)),
        (2, [3, 5]),
        (3, (2, 3, 4)),
    )

    @testset "N=$N, max_occupancies=$max_occupancies" for (N, max_occupancies) in test_cases
        max_fock_states = maximum(max_occupancies) + 1

        fock_space    = FockSpace(Val(N), max_occupancies)
        bcf           = random_bcf(N)
        atom_params   = AtomParams(rand(), rand(ComplexF64), rand(ComplexF64))
        solver_params = create_solver_params(bcf, fock_space, atom_params)

        hme! = HME{N,max_fock_states}()

        @test hme! isa HME{N,max_fock_states}

        @testset "Type parameters" begin
            @test typeof(hme!).parameters[1] == N
            @test typeof(hme!).parameters[2] == max_fock_states
        end

        @testset "Fields" begin
            @test length(hme!._sqrt_of_int) == max_fock_states

            for idx in eachindex(hme!._sqrt_of_int)
                @test hme!._sqrt_of_int[idx] == sqrt(idx - 1)
            end

            @test hme!.f_tmp isa MVector{N,ComplexF64}
            @test hme!.g_tmp isa MVector{N,ComplexF64}
        end
    end

    @testset "Invalid inputs" begin
        @test_throws AssertionError HME{0,1}()
        @test_throws AssertionError HME{1,0}()
        @test_throws AssertionError HME{1,-1}()
        @test_throws AssertionError HME{-1,1}()
        @test_throws AssertionError HME{1.0,1}()
        @test_throws AssertionError HME{1,1.0}()
    end
end;

@testset "HME (function)" begin
    println("HME: evaluation")
    test_cases = (
        (1, 0),
        (1, (5,)),
        (2, [3, 5]),
        (3, (2, 3, 4)),
    )

    @testset "N=$N, max_occupancies=$max_occupancies" for (N, max_occupancies) in test_cases
        max_fock_states = maximum(max_occupancies) + 1

        fock_space = FockSpace(Val(N), max_occupancies)
        bcf = random_bcf(N)
        atom_params = AtomParams(
            rand(),
            rand(ComplexF64),
            rand(ComplexF64),
        )
        solver_params = create_solver_params(
            bcf,
            fock_space,
            atom_params,
        )

        hme! = HME{N,max_fock_states}()

        fock_dim = solver_params.fock_dim
        t = rand()

        ρ = rand(ComplexF64, 2, 2, fock_dim, fock_dim)
        dρ = rand(ComplexF64, 2, 2, fock_dim, fock_dim)

        @testset "Type stability" begin
            @inferred hme!(dρ, ρ, solver_params, t)
        end

        @testset "Evaluation" begin
            hme!(dρ, ρ, solver_params, t)

            @test all(isfinite, real.(dρ))
            @test all(isfinite, imag.(dρ))
        end

        @testset "Mutation" begin
            @. dρ = 0

            ρ_copy = deepcopy(ρ)
            dρ_copy = deepcopy(dρ)
            solver_params_copy = deepcopy(solver_params)

            hme!(dρ_copy, ρ_copy, solver_params_copy, t)

            @test size(dρ_copy) == (2, 2, fock_dim, fock_dim)

            @test ρ_copy == ρ
            @test solver_params_copy == solver_params
            @test dρ_copy != dρ

            hme!(dρ, ρ_copy, solver_params_copy, t)

            @test dρ == dρ_copy
        end
    end

    @testset "Allocations" begin
        N = 3
        max_occupancies = (0,0,0)

        max_fock_states = maximum(max_occupancies) + 1

        fock_space    = FockSpace(Val(N), max_occupancies)
        bcf           = random_bcf(N)
        atom_params   = AtomParams(rand(), rand(ComplexF64), rand(ComplexF64))
        solver_params = create_solver_params(bcf, fock_space, atom_params)

        hme! = HME{N,max_fock_states}()

        fock_dim = solver_params.fock_dim
        t = rand()

        ρ = rand(ComplexF64, 2, 2, fock_dim, fock_dim)
        dρ = similar(ρ)

        @test @ballocated($hme!($dρ, $ρ, $solver_params, $t)) == 0
    end

    println("HME: done")
end;