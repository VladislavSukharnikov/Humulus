# =============================================================================
# Solver parameters. 
# =============================================================================

"""
    create_solver_params(bcf, fock_space, atom_params)
Construct a `NamedTuple` containing the parameters required to evaluate the
right-hand side of the equations of motion.

The returned object is intended for repeated calls to `solve_hme` and
`solve_hops`.

# Arguments
- `bcf::BCF`: bath correlation function.
- `fock_space::FockSpace`: pseudo-Fock space.
- `atom_params::AtomParams`: atomic parameters.

# Returns
A `NamedTuple` containing:
- atomic parameters;
- bath-correlation-function parameters;
- pseudo-Fock-space data.
"""
function create_solver_params(
                            bcf::BCF, 
                            fock_space::FockSpace, 
                            atom_params::AtomParams
                          )::NamedTuple
                              
    @assert n_modes(fock_space)==n_modes(bcf) "Incompatible number of modes: FockSpace has $(n_modes(fock_space)), BCF has $(n_modes(bcf))."
    
    N = n_modes(bcf)

    (; ν₀, c_g, c_e) = atom_params
    (; Γ, G, f_vector, g_vector) = bcf
    (; basis_states, fock_dim, max_occupancies, raise_index, lower_index) = fock_space

    return (; N, 
              ν₀, c_g, c_e, 
              Γ, G, f_vector, g_vector,
              basis_states, fock_dim, max_occupancies, raise_index, lower_index)
end