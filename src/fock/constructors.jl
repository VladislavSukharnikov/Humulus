"""
    FockSpace(::Val{N}, max_occupancies; KeyType=Int, IndType=Int)

Construct and return a truncated `N`-mode Fock space.

The truncation is specified by `max_occupancies`, which may be

- an `Int`, retaining all basis states satisfying `∑ᵢ nᵢ ≤ max_occupancies`,
- an `NTuple{N,Int}`, specifying an independent occupation cutoff for each mode, or
- a `Vector{Int}` of length `N`, which is converted internally to an `NTuple`.

The constructor generates the Fock-space basis together with lookup tables
for basis indexing and neighboring basis states.

# Arguments
- `max_occupancies`: occupation-number truncation.

# Keyword arguments
- `KeyType`: integer type used to store occupation numbers.
- `IndType`: integer type used to store basis indices.

# Exceptions

An exception is thrown if

- `N ≤ 0`;
- any truncation level is negative;
- the number of truncation levels does not equal `N`; or
- the resulting Fock space cannot be represented using the selected `KeyType` or `IndType`.
"""
function FockSpace(
                ::Val{N},
                max_occupancies::Union{NTuple{N,Int},Int,Vector{Int}};
                KeyType::Type{<:Integer}=Int,
                IndType::Type{<:Integer}=Int
            ) where {N}
        
    # Input validation.    
    N > 0 || throw(DomainError(N, "The number of modes must be positive."))
    if max_occupancies isa Union{Tuple, Vector{Int}}

        length(max_occupancies) == N ||
            throw(DimensionMismatch("Expected $N truncation levels, but got $(length(max_occupancies))."))
        all(>=(0), max_occupancies) ||
            throw(DomainError(max_occupancies, "All truncation levels must be non-negative."))

        fock_dim = prod(big.(max_occupancies .+ 1))
        fock_dim <= typemax(Int) ||
            throw(ArgumentError("The resulting Fock space contains $fock_dim states, which exceeds the maximum array index size ($(typemax(Int)))."))
        fock_dim <= typemax(IndType) ||
            throw(ArgumentError("The selected IndType cannot represent a Fock space of size $fock_dim."))

        all(<=(typemax(KeyType)), max_occupancies) ||
            throw(ArgumentError("The selected KeyType cannot represent one or more truncation levels."))

    elseif max_occupancies isa Integer

        max_occupancies >= 0 ||
            throw(DomainError(max_occupancies, "The truncation level must be non-negative."))
        max_occupancies <= typemax(KeyType) ||
            throw(ArgumentError("The selected KeyType cannot represent the truncation level $max_occupancies."))

        fock_dim = binomial(big(max_occupancies + N), N)
        fock_dim <= typemax(Int) ||
            throw(ArgumentError("The resulting Fock space contains $fock_dim states, which exceeds the maximum array index size ($(typemax(Int)))."))
        fock_dim <= typemax(IndType) ||
            throw(ArgumentError("The selected IndType cannot represent a Fock space of size $fock_dim."))

    else
        throw(ArgumentError("`max_occupancies` must be either an Integer or a tuple/vector of integers."))
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
    @assert length(state_to_index) == length(basis_states) == fock_dim  "Internal error: inconsistent Fock space size."
    
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