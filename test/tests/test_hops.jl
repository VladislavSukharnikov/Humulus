@testset "HOPS" begin
    @info "HOPS: construction..."

    test_cases = (
        (1, 0),
        (1, (5,)),
        (2, (3, 5)),
        (3, (2, 3, 4)),
    )

    time_grid = TimeGrid(0.0, 20.0, 200)

    @testset "N=$N, max_occupancies=$max_occupancies" for
        (N, max_occupancies) in test_cases

        bcf = random_bcf(N)

        max_fock_states = maximum(max_occupancies) + 1
        hops! = HOPS{N,max_fock_states}(time_grid)

        @test hops! isa HOPS{N,max_fock_states}

        @testset "Type parameters" begin
            @test typeof(hops!).parameters[1] == N
            @test typeof(hops!).parameters[2] == max_fock_states
        end

        @testset "Fields" begin
            @test length(hops!._sqrt_of_int) == max_fock_states

            for idx in eachindex(hops!._sqrt_of_int)
                @test hops!._sqrt_of_int[idx] == sqrt(idx - 1)
            end

            @test hops!.f_tmp isa MVector{N,ComplexF64}
            @test hops!.g_tmp isa MVector{N,ComplexF64}
            @test hops!._time_grid == time_grid
            @test hops!.noise isa Vector{ComplexF64}
            @test length(hops!.noise) == length(hops!._time_grid)
        end
    end

    @testset "Invalid inputs" begin

        @test_throws ArgumentError HOPS{1.0,1}(time_grid)
        @test_throws ArgumentError HOPS{1,1.0}(time_grid)

        @test_throws DomainError HOPS{0,1}(time_grid)
        @test_throws DomainError HOPS{1,0}(time_grid)
        @test_throws DomainError HOPS{1,-1}(time_grid)
        @test_throws DomainError HOPS{-1,1}(time_grid)
    end
end;


@testset "HOPS (function)" begin
    @info "HOPS: evaluation..."

    test_cases = (
        (1, 0),
        (1, (5,)),
        (2, (3, 5)),
        (3, (2, 3, 4)),
    )

    time_grid = TimeGrid(0.0, 20.0, 200)

    @testset "N=$N, max_occupancies=$max_occupancies" for
        (N, max_occupancies) in test_cases

        bcf = random_bcf(N)

        atom_params     = AtomParams(rand(), rand(ComplexF64), rand(ComplexF64))
        fock_space      = FockSpace(Val(N), max_occupancies)
        max_fock_states = maximum(max_occupancies) + 1
        hops! = HOPS{N,max_fock_states}(time_grid)

        @test all(iszero, hops!.noise)
        for i in eachindex(hops!.noise)
            hops!.noise[i] = rand(ComplexF64)
        end
        @test all(!iszero, hops!.noise)

        fock_dim = fock_space.fock_dim

        dψ = rand(ComplexF64, 2, fock_dim)
        ψ  = rand(ComplexF64, 2, fock_dim)
        dX = rand(ComplexF64, N, 1)
        X  = rand(ComplexF64, N, 1)
        du = ArrayPartition(dψ, dX)
        u  = ArrayPartition(ψ, X)

        solver_params = create_solver_params(bcf, fock_space, atom_params)

        t = rand() * time_grid.t_end

        @testset "Type stability" begin
            @inferred hops!(du, u, solver_params, t)
        end

        @testset "Evaluation" begin
            hops!(du, u, solver_params, t)

            @test all(isfinite, real.(du.x[1]))
            @test all(isfinite, imag.(du.x[1]))
            @test all(isfinite, real.(du.x[2]))
            @test all(isfinite, imag.(du.x[2]))
        end

        @testset "Mutation" begin
            @. du = 0.0

            u_copy  = deepcopy(u)
            du_copy = deepcopy(du)
            solver_params_copy = deepcopy(solver_params)

            hops!(du_copy, u_copy, solver_params_copy, t)

            @test u_copy == u
            @test solver_params_copy == solver_params
            @test du_copy != du

            hops!(du, u_copy, solver_params_copy, t)

            @test du == du_copy
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

        hops! = HOPS{N,max_fock_states}(time_grid)

        fock_dim = fock_space.fock_dim

        dψ = rand(ComplexF64, 2, fock_dim)
        ψ  = rand(ComplexF64, 2, fock_dim)
        dX = rand(ComplexF64, N, 1)
        X  = rand(ComplexF64, N, 1)
        du = ArrayPartition(dψ, dX)
        u  = ArrayPartition(ψ, X)

        solver_params = create_solver_params(bcf, fock_space, atom_params)

        t = rand() * time_grid.t_end

        for i in eachindex(hops!.noise)
            hops!.noise[i] = rand(ComplexF64)
        end
        
        hops!(du, u, solver_params, t)

        @test @ballocated($hops!($du, $u, $solver_params, $t)) == 0
    end

    @info "HOPS: completed."
end;