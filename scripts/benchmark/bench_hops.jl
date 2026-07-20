include("bench_utils.jl");

# =============================================================================
# Benchmark HOPS construction.
# =============================================================================

let
    title = "Benchmark HOPS construction."
    print_header(title)

    N = 5

    # Maximum occupation number for each mode.
    max_occupancies = 100

    # Maximum number of Fock states retained across all modes.
    max_fock_states = maximum(max_occupancies) + 1

    # Noise discretization grid.
    time_grid = Humulus.TimeGrid(0.0, 1.0, 10_000)

    @info "Benchmark parameters" N max_occupancies time_grid.n_points

    # Construct and benchmark HOPS.
    hops! = Humulus.HOPS{N,max_fock_states}(time_grid)
    
    bench = @benchmark Humulus.HOPS{$N,$max_fock_states}($time_grid)
    benchmark_construction(hops!, bench)
end


# =============================================================================
# Benchmark HME evaluation.
# =============================================================================

let
    title = "Benchmark HOPS evaluation."
    print_header(title)

    N = 2

    # Construct BCF.
    bcf = Humulus.random_bcf(N)

    # Construct AtomParams.
    atom_params = AtomParams(rand(), rand(ComplexF64), rand(ComplexF64))

    # Construct FockSpace.
    max_occupancies = 50
    fock_space      = Humulus.FockSpace(Val(N), max_occupancies, Int, Int)
    max_fock_states = maximum(max_occupancies) + 1

    # Noise discretization grid.
    time_grid = Humulus.TimeGrid(0.0, 1.0, 10_000)

    # Construct HOPS.
    hops! = Humulus.HOPS{N,max_fock_states}(time_grid)

    # Generate arguments for hops!.
    fock_dim = fock_space.fock_dim
    dψ = rand(ComplexF64, 2, fock_dim)
    ψ  = rand(ComplexF64, 2, fock_dim)
    dX = rand(ComplexF64, N, 1)
    X  = rand(ComplexF64, N, 1)
    du = ArrayPartition(dψ, dX)
    u  = ArrayPartition(ψ, X)
    solver_params = Humulus.create_solver_params(bcf, fock_space, atom_params)
    t = rand()

    # Trigger compilation.
    hops!(du, u, solver_params, t)

    @info "Benchmark parameters:" N fock_dim time_grid.n_points

    bench = @benchmark $hops!($du, $u, $solver_params, $t)
    benchmark_evaluation(bench)
end