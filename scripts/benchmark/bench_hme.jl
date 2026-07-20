include("bench_utils.jl");

# =============================================================================
# Benchmark HME construction.
# =============================================================================

let
    title = "Benchmark HME construction."
    print_header(title)

    N = 5

    # Maximum occupation number for each mode.
    max_occupancies = [10, 10, 10, 10, 10]

    # Maximum number of Fock states retained across all modes.
    max_fock_states = maximum(max_occupancies) + 1

    @info "Benchmark parameters" N max_occupancies

    # Construct and benchmark HME.
    hme!  = Humulus.HME{N,max_fock_states}()
    
    bench = @benchmark Humulus.HME{$N,$max_fock_states}()
    benchmark_construction(hme!, bench)
end


# =============================================================================
# Benchmark HME evaluation.
# =============================================================================

let
    title = "Benchmark HME evaluation."
    print_header(title)

    N = 2

    # Construct BCF.
    bcf = Humulus.random_bcf(N)

    # Construct AtomParams.
    atom_params = AtomParams(rand(), rand(ComplexF64), rand(ComplexF64))

    # Construct FockSpace.
    max_occupancies = (10, 10)
    fock_space      = Humulus.FockSpace(Val(N), max_occupancies, Int, Int)
    max_fock_states = maximum(max_occupancies) + 1

    # Construct HME.
    hme! = Humulus.HME{N,max_fock_states}()

    # Generate arguments for hme!.
    fock_dim = fock_space.fock_dim
    dρ = rand(ComplexF64, 2, 2, fock_dim, fock_dim)
    ρ  = rand(ComplexF64, 2, 2, fock_dim, fock_dim)
    solver_params = Humulus.create_solver_params(bcf, fock_space, atom_params)
    t = rand()

    # Trigger compilation.
    hme!(dρ, ρ, solver_params, t)

    @info "Benchmark parameters:" N fock_dim

    bench = @benchmark $hme!($dρ, $ρ, $solver_params, $t)
    benchmark_evaluation(bench)
end