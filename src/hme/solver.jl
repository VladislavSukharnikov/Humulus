# =============================================================================
# Time evolution
# =============================================================================

"""
    solve_hme(grid_params, bcf, atom_params, max_occupancies;
              KeyType=Int, IndType=Int)

Solve the hierarchy of master equations (HME).

This method constructs the pseudo-Fock space and solver parameters before
integrating the HME.

# Arguments
- `grid_params::GridParams`: time discretization parameters.
- `bcf::BCF`: bath correlation function.
- `atom_params::AtomParams`: atomic parameters and initial state.
- `max_occupancies::Union{NTuple{N,Int},Int,Vector{Int}}`: pseudo-Fock-space truncation.

# Keyword arguments
- `KeyType::Type{<:Integer}`: integer type used for pseudo-Fock-space keys.
- `IndType::Type{<:Integer}`: integer type used for pseudo-Fock-space indices.

# Exceptions

Any exceptions thrown while constructing the pseudo-Fock space or solver
parameters are propagated.
"""
function solve_hme(
                grid_params::GridParams,
                bcf::BCF{N},
                atom_params::AtomParams,
                max_occupancies::Union{NTuple{N,Int},Int,Vector{Int}};
                KeyType::Type{<:Integer}=Int,
                IndType::Type{<:Integer}=Int,
            ) where {N}
    
    # Construct the pseudo-Fock space and solver parameters.
    fock_space      = FockSpace(Val(N), max_occupancies, KeyType, IndType)
    solver_params   = create_solver_params(bcf, fock_space, atom_params)
    max_fock_states = maximum(max_occupancies)+1

    # Dispatch to the implementation specialized on the problem size.
    return _solve_hme(
                    Val(N), 
                    Val(max_fock_states), 
                    grid_params,
                    solver_params,
                )
end


"""
    solve_hme(grid_params, solver_params)

Solve the hierarchy of master equations (HME).

The `solver_params` argument should be created with
`create_solver_params`.

# Arguments
- `grid_params::GridParams`: time discretization parameters.
- `solver_params`: solver parameters returned by `create_solver_params`.
"""
function solve_hme(
                grid_params::GridParams,
                solver_params::NamedTuple,
            )

    (; N, max_occupancies) = solver_params
    max_fock_states = maximum(max_occupancies) + 1

    # Dispatch to the implementation specialized on the problem size.
    return _solve_hme(
                    Val(N),
                    Val(max_fock_states),
                    grid_params,
                    solver_params,
                )
end


# Internal implementation specialized on the number of modes and the maximum
# number of retained Fock states.
@inline function _solve_hme(
                    ::Val{N},
                    ::Val{MaxFockStates},
                    grid_params::GridParams,
                    solver_params::NamedTuple, 
                ) where {N,MaxFockStates}
    
    (; ts_save, dt_max) = grid_params;

    # Allocate the output density matrix.
    ρ_s = zeros(ComplexF64, 2, 2, ts_save.n_points)

    # Construct the HME right-hand-side functor and ODE problem.
    (; fock_dim, c_g, c_e) = solver_params
    hme! = HME{N,MaxFockStates}();
    ρ_0  = init_hme(fock_dim, c_g, c_e)
    prob = ODEProblem{true}(ODEFunction{true}(hme!), ρ_0, (ts_save.t_start, ts_save.t_end), solver_params);

    # Solve ODE.
    data = solve(
            prob, 
            Tsit5(), 
            adaptive=true, 
            dtmax=dt_max, 
            saveat=ts_save, 
            dense = false, 
            save_everystep = false, 
            calck = false,
            progress = true,
        )

    # Save data.
    for t_idx in 1:ts_save.n_points
        ρ_s_slice = @view ρ_s[:,:,t_idx]
        u_slice = @view data[:, :, 1, 1, t_idx]
        ρ_s_slice .= u_slice
    end

    return ρ_s
end