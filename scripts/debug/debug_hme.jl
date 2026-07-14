include("debug_utils.jl");

# =============================================================================
# Debugging: @code_warntype inspection of hme!(...).
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

    # Construct HME.
    hme! = Humulus.HME{N,max_fock_states}();

    # Generate arguments for hme!.
    fock_dim = fock_space.fock_dim
    dρ = rand(ComplexF64, 2, 2, fock_dim, fock_dim)
    ρ  = rand(ComplexF64, 2, 2, fock_dim, fock_dim)
    solver_params = Humulus.create_solver_params(bcf, fock_space, atom_params)
    t = rand()

    # Trigger compilation.
    hme!(dρ, ρ, solver_params, t)

    @code_warntype hme!(dρ, ρ, solver_params, t)
end


# =============================================================================
# Debugging: @code_warntype inspection of _solve_hme(...).
# =============================================================================

let 
    # Construct BCF.
    N = 2;
    bcf = Humulus.random_bcf(N)
    
    # Construct AtomParams.
    atom_params   = AtomParams(rand(), rand(ComplexF64), rand(ComplexF64))
    
    # Construct GridParams.
    grid_params = GridParams(1.0, 100, 1)

    # Construct FockSpace.
    max_occupancies  = (10, 10);
    KeyType, IndType = Int16, Int32;
    fock_space       = Humulus.FockSpace(Val(N), max_occupancies, KeyType, IndType)
    max_fock_states  = maximum(max_occupancies) + 1

    # Create solver parameters.
    solver_params = Humulus.create_solver_params(bcf, fock_space, atom_params)

    # Trigger compilation.
    Humulus._solve_hme(
                        Val(N),
                        Val(max_fock_states),
                        grid_params,
                        solver_params, 
                    )
                    
    @code_warntype Humulus._solve_hme(
                        Val(N),
                        Val(max_fock_states),
                        grid_params,
                        solver_params, 
                    )
end