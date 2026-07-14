# =============================================================================
# Time evolution
# =============================================================================

"""
    solve_hme(grid_params, bcf, atom_params, max_occupancies; KeyType=Int, IndType=Int)

Solve the hierarchy of master equations (HME).

This method constructs the pseudo-Fock space and solver parameters before
integrating the HME.

# Arguments
- `grid_params::GridParams`: time discretization parameters.
- `bcf::BCF`: bath correlation function.
- `atom_params::AtomParams`: atomic parameters and initial state.
- `max_occupancies::Union{NTuple{N,Int},Int}`: pseudo-Fock-space truncation.

# Keyword arguments
- `KeyType::Type{<:Integer}`: integer type used for pseudo-Fock-space keys.
- `IndType::Type{<:Integer}`: integer type used for pseudo-Fock-space indices.

# Returns
- `Array{ComplexF64,3}`: reduced density matrix computed on the save-time grid.
"""
function solve_hme(
                grid_params::GridParams, 
                bcf::BCF{N}, 
                atom_params::AtomParams, 
                max_occupancies::Union{NTuple{N,Int},Int};
                KeyType::Type{<:Integer}=Int,
                IndType::Type{<:Integer}=Int,
            ) where {N}

    fock_space      = FockSpace(Val(N), max_occupancies, KeyType, IndType)
    solver_params   = create_solver_params(bcf, fock_space, atom_params)
    max_fock_states = maximum(max_occupancies)+1

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

# Returns
- `Array{ComplexF64,3}`: reduced density matrix computed on the save-time grid.
"""
function solve_hme(
                grid_params::GridParams,
                solver_params::NamedTuple,
            )

    (; N, max_occupancies) = solver_params
    max_fock_states = maximum(max_occupancies) + 1

    return _solve_hme(
                    Val(N),
                    Val(max_fock_states),
                    grid_params,
                    solver_params,
                )
end


# Internal realization.
@inline function _solve_hme(
                    ::Val{N},
                    ::Val{MaxFockStates},
                    grid_params::GridParams,
                    solver_params::NamedTuple, 
                ) where {N,MaxFockStates}
    
    (; ts_save, dt_max) = grid_params;

    # Container for the physical density matrix
    ρ_s = zeros(ComplexF64, 2, 2, ts_save.n_points)

    (; fock_dim, c_g, c_e) = solver_params
    hme! = HME{N,MaxFockStates}();
    ρ_0  = init_hme(fock_dim, c_g, c_e)
    prob = ODEProblem{true}(ODEFunction{true}(hme!), ρ_0, (ts_save.t_start, ts_save.t_end), solver_params);

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

    for save_idx in 1:ts_save.n_points
        ρ_s_slice = @view ρ_s[:,:,save_idx]
        u_slice = view(data.u[save_idx],:,:,1,1)
        ρ_s_slice .= u_slice
    end

    return ρ_s
end
