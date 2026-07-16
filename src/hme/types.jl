# =============================================================================
# Hierarchy of Master Equations (HME) type.
# =============================================================================

"""
    HME{N,MaxFockStates} <: Function

Functor implementing the hierarchy of master equations (HME).

The object stores reusable temporary buffers and a lookup table of square
roots used when applying creation and annihilation operators.

# Type parameters
- `N`: number of effective modes.
- `MaxFockStates`: maximum number of retained Fock states per mode.

# Fields
- `f_tmp`: temporary storage for the values of the `f` functions.
- `g_tmp`: temporary storage for the values of the `g` functions.
- `_sqrt_of_int`: lookup table containing `sqrt.(0:MaxFockStates-1)`.
"""
mutable struct HME{N,MaxFockStates}<:Function
    f_tmp::MVector{N,ComplexF64}
    g_tmp::MVector{N,ComplexF64}
    _sqrt_of_int::SVector{MaxFockStates,Float64}   
end


"""
    HME{N,MaxFockStates}()

Construct an `HME` functor for an `N`-mode system.

`MaxFockStates` specifies the maximum number of retained Fock states per mode.
"""
function HME{N,MaxFockStates}() where {N,MaxFockStates}

    # Input validation. 
    N isa Int ||
        throw(ArgumentError("The number of modes `N` must be an `Int`, got $(typeof(N))."))

    MaxFockStates isa Int ||
        throw(ArgumentError("`MaxFockStates` must be an `Int`, got $(typeof(MaxFockStates))."))

    N > 0 ||
        throw(DomainError(N, "The number of modes `N` must be positive."))

    MaxFockStates >= 1 ||
        throw(DomainError(MaxFockStates, "`MaxFockStates` must be at least 1, since the vacuum state is always included."))

    # Allocate temporary work buffers.
    f_tmp = zeros(MVector{N,ComplexF64})
    g_tmp = zeros(MVector{N,ComplexF64})

    # Precompute square roots of the occupation numbers.
    sqrt_of_int = SVector{MaxFockStates,Float64}(sqrt(i) for i in 0:MaxFockStates-1)

    # Construct the HME object.
    return HME{N,MaxFockStates}(f_tmp, g_tmp, sqrt_of_int)
end