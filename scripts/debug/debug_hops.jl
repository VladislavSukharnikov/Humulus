include("debug_utils.jl");

# =============================================================================
# Debugging: @code_warntype inspection of hops!(...).
# =============================================================================

let 
    # Construct BCF.
    N = 2
    bcf = Humulus.random_bcf(N)
    
    # Construct AtomParams.
    atom_params = AtomParams(rand(), rand(ComplexF64), rand(ComplexF64))
    
    # Construct FockSpace.
    max_occupancies  = (10, 10);
    KeyType, IndType = Int16, Int32;
    fock_space       = Humulus.FockSpace(Val(N), max_occupancies, KeyType, IndType)
    max_fock_states  = maximum(max_occupancies) + 1

    # Construct HOPS.
    time_grid = Humulus.TimeGrid(0.0, 1.0, 100)
    hops! = Humulus.HOPS{N,max_fock_states}(time_grid)

    # Generate arguments for hme!.
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

    @code_warntype hops!(du, u, solver_params, t)
end


# =============================================================================
# Debugging: @code_warntype inspection of _solve_hops(...).
# =============================================================================

let 
    # Construct BCF.
    N = 2
    bcf = Humulus.random_bcf(N)
    
    # Construct AtomParams.
    atom_params = AtomParams(rand(), rand(ComplexF64), rand(ComplexF64))
    
    # Construct FockSpace.
    max_occupancies  = (10, 10);
    KeyType, IndType = Int16, Int32;
    fock_space       = Humulus.FockSpace(Val(N), max_occupancies, KeyType, IndType)
    max_fock_states  = maximum(max_occupancies) + 1

    # Construct GridParams.
    grid_params = GridParams(1.0, 100, 1)

    # Create solver parameters.
    solver_params = Humulus.create_solver_params(bcf, fock_space, atom_params)

    # Create noise cache.
    (; dt_max, ts_save) = grid_params
    noise_oversampling = 2
    dt_noise = dt_max / noise_oversampling
    grid_size = ceil(Int, (ts_save.t_end) / dt_noise)
    path = Humulus.get_bcf_eigen_cache(bcf, ts_save.t_end, grid_size)

    n_trajectories = 1

    # Trigger compilation.
    output = Humulus._solve_hops(
                    Val(N), 
                    Val(max_fock_states),  
                    grid_params, 
                    solver_params, 
                    n_trajectories,
                    path,
                )

    @code_warntype Humulus._solve_hops(
                    Val(N), 
                    Val(max_fock_states),  
                    grid_params, 
                    solver_params, 
                    n_trajectories,
                    path,
                )
    rm(path)
end