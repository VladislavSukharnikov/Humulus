include("prof_utils.jl")

# =============================================================================
# Profile repeated BCF evaluation for the three-mode squeezed reservoir.
# =============================================================================

let
    # Number of modes.
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
    profile_hme_eval!(hme!, dρ, ρ, solver_params, t);

    repetition_number = 100;
    @profview profile_hme_eval!(hme!, dρ, ρ, solver_params, t; repetition_number=repetition_number)
end
