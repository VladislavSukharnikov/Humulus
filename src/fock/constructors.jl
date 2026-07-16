"""
    FockSpace(::Val{N}, max_occupancies; KeyType=Int, IndType=Int)

Construct a truncated `N`-mode Fock space.

The truncation is specified by `max_occupancies`, which may be either

- an `Int`, retaining all basis states satisfying `∑ᵢ nᵢ ≤ max_occupancies`, or
- an `NTuple{N,Int}`, specifying an independent occupation cutoff for each mode.

The constructor generates the Fock-space basis together with lookup tables
for basis indexing and neighboring basis states.

# Arguments
- `max_occupancies`: occupation-number truncation.

# Keyword arguments
- `KeyType`: integer type used to store occupation numbers.
- `IndType`: integer type used to store basis indices.

# Returns
- `FockSpace`: truncated Fock-space representation with precomputed lookup tables.
"""
function FockSpace(
                ::Val{N},
                max_occupancies::Union{NTuple{N,Int},Int,Vector{Int}},
                KeyType::Type{<:Integer}=Int,
                IndType::Type{<:Integer}=Int
            ) where {N}

    @assert N isa Int "Number of modes must be integer."
    @assert N>0       "Number of modes must be positive."

    if max_occupancies isa Tuple || max_occupancies isa Vector{Int} 
        fock_dim = prod(big.(max_occupancies.+1))
        @assert length(max_occupancies) == N "Incorrect size of the truncation list."
        @assert fock_dim <= typemax(Int) "The Fock space size exceeds the maximum array index size."
        @assert all(>=(0), max_occupancies) "All truncation levels must be non-negative."
        @assert fock_dim <= typemax(IndType) "The Fock space size exceeds the range of `IndType`."
        @assert all(<=(typemax(KeyType)), max_occupancies)  "Some truncation levels exceed the range of `KeyType`."
    elseif max_occupancies isa Integer
        fock_dim = binomial(big(max_occupancies+N),N)
        @assert fock_dim <= typemax(Int) "The Fock space size exceeds the maximum array index size."
        @assert max_occupancies >= 0 "The truncation level must be non-negative."
        @assert fock_dim <= typemax(IndType) "The Fock space size exceeds the range of `IndType`."
        @assert max_occupancies<=typemax(KeyType) "Some truncation levels exceed the range of `KeyType`."
    end

    if max_occupancies isa Vector
        max_occupancies = NTuple{N,Int}(max_occupancies)
    end
    
    fock_dim  = Int(fock_dim)

    # Allocate storage for the basis lookup tables.
    state_to_index = Dict{SVector{N,KeyType}, IndType}()
    sizehint!(state_to_index, fock_dim)

    # Allocate storage all basis states.
    basis_states = SVector{N,KeyType}[]
    sizehint!(basis_states, fock_dim)
    
    # Create mapping from occupation number sets to integers.
    build_basis!(
        state_to_index, 
        basis_states, 
        Val(N), 
        max_occupancies, 
        KeyType, 
        IndType,
    )

    # Verify that the generated basis has the expected size.
    @assert length(state_to_index)==length(basis_states)==fock_dim "The generated Fock space has an unexpected size."

    # Precompute neighboring basis indices for raising and lowering operators.
    raise_index = zeros(IndType, N, fock_dim)
    lower_index = zeros(IndType, N, fock_dim)

    @inbounds for index in 1:fock_dim
        state = basis_states[index]
        for j in 1:N
            lower_index[j,index] = neighbor_index(state_to_index, state, j, -1)
            raise_index[j,index] = neighbor_index(state_to_index, state, j, +1)
        end
    end

    return FockSpace{N,KeyType,IndType}(
                max_occupancies, 
                state_to_index, 
                basis_states, 
                fock_dim, 
                raise_index, 
                lower_index,
            )
end