# =============================================================================
# Basis generation.
# =============================================================================

"""
    build_basis!(
        state_to_index,
        basis_states,
        Val(N),
        max_occupancies,
        KeyType,
        IndType,
    )

Populate `basis_states` and `state_to_index` with the basis of an `N`-mode
Fock space truncated by per-mode occupation cutoffs.

Each basis state is an `SVector{N,KeyType}` satisfying
`0 ≤ nᵢ ≤ max_occupancies[i]` for `i = 1, …, N`.

The basis is generated in lexicographic order using `Iterators.product`.
`basis_states` stores the basis states, while `state_to_index` maps each
basis state to its corresponding 1-based index.

This function mutates `state_to_index` and `basis_states` and returns
`nothing`.
"""
function build_basis!(
                    state_to_index::Dict{SVector{N,KeyType},IndType},
                    basis_states::Vector{SVector{N,KeyType}},
                    ::Val{N},
                    max_occupancies::NTuple{N,Int}, 
                    ::Type{KeyType}, 
                    ::Type{IndType}
                ) where {N,KeyType<:Integer,IndType<:Integer}

    ranges = SVector{N, UnitRange{Int64}}(KeyType(0):KeyType(max_val) for max_val in max_occupancies)
    index  = IndType(1)
    for combination in Iterators.product(ranges...)
        state_to_index[SVector{N,KeyType}(combination)] = index
        push!(basis_states, combination)
        index = IndType(index+1)
    end
    return nothing
end


"""
    build_basis!(
        state_to_index,
        basis_states,
        Val(N),
        max_occupancies,
        KeyType,
        IndType,
    )

Populate `basis_states` and `state_to_index` with the basis of an `N`-mode
Fock space truncated by a global occupation-number cutoff.

Only basis states satisfying

    ∑ᵢ nᵢ ≤ max_occupancies

are included.

The basis is generated in lexicographic order using `Iterators.product`.
`basis_states` stores the basis states, while `state_to_index` maps each
basis state to its corresponding 1-based index.

This function mutates `state_to_index` and `basis_states` and returns
`nothing`.
"""
function build_basis!(
                    state_to_index::Dict{SVector{N,KeyType},IndType},
                    basis_states::Vector{SVector{N,KeyType}},
                    ::Val{N},
                    max_occupancies::Int, 
                    ::Type{KeyType}, 
                    ::Type{IndType}
                ) where {N,KeyType<:Integer,IndType<:Integer}

    ranges = SVector{N,UnitRange{Int64}}(KeyType(0):KeyType(max_occupancies) for _ in 1:N)
    index  = IndType(1)
    for combination in Iterators.product(ranges...)
        if sum(combination) <= max_occupancies
            state_to_index[SVector{N, KeyType}(combination)] = index
            push!(basis_states, combination)
            index = IndType(index+1)
        end
    end
    return nothing
end


# =============================================================================
# Basis lookup. 
# =============================================================================

"""
    neighbor_index(state_to_index, state, mode, occupation_shift)

Return the basis index of the state obtained by shifting the occupation
number of `state` in `mode` by `occupation_shift`.

# Arguments
- `state_to_index`: mapping from basis states to basis indices.
- `state`: occupation-number state.
- `mode`: mode whose occupation number is modified.
- `occupation_shift`: change in the occupation number.

# Returns

Returns the corresponding basis index, or zero if the resulting state is
not contained in the truncated Fock-space basis.
"""
@inline function neighbor_index(
                        state_to_index::Dict{SVector{N,KeyType},IndType}, 
                        state::SVector{N,KeyType}, 
                        mode::Int, 
                        occupation_shift::Integer
                    ) where {N,KeyType<:Integer,IndType<:Integer}

    # Reject shifts that produce invalid occupation numbers.
    @assert 1 <= mode <= N "Internal error: invalid mode index $mode (valid range is 1:$N)."
    
    if occupation_shift<0 && state[mode]<abs(occupation_shift)
        return zero(IndType)
    elseif state[mode]+occupation_shift > typemax(KeyType)
        return zero(IndType)
    end

    # Look up the neighboring basis state.
    neighbor_state = setindex(state, KeyType(state[mode]+occupation_shift), mode)

    return get(state_to_index, neighbor_state, zero(IndType))
end