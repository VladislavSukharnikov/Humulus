# =============================================================================
# Fock space type. 
# =============================================================================

"""
    FockSpace{N,KeyType,IndType}

Truncated Fock space for an `N`-mode system.

The type stores the basis states together with lookup tables for converting
between occupation-number states and basis indices, as well as precomputed
neighbor indices for creation and annihilation operators.

# Type parameters
- `N`: number of modes.
- `KeyType`: integer type used to store occupation numbers.
- `IndType`: integer type used to store basis indices.

# Fields
- `max_occupancies`: maximum occupation number for each mode.
- `state_to_index`: mapping from occupation-number states to basis indices.
- `basis_states`: basis states ordered by basis index.
- `fock_dim`: dimension of the truncated Fock space.
- `raise_index`: basis indices of states reached by applying a creation operator.
- `lower_index`: basis indices of states reached by applying an annihilation operator.
"""
struct FockSpace{N,KeyType<:Integer,IndType<:Integer}
    max_occupancies :: Union{Int,NTuple{N,Int}}
    state_to_index  :: Dict{SVector{N,KeyType},IndType}
    basis_states    :: Vector{SVector{N,KeyType}}
    fock_dim        :: Int
    raise_index     :: Matrix{IndType}
    lower_index     :: Matrix{IndType}
end

# Number of effective modes.
n_modes(::FockSpace{N}) where {N} = N