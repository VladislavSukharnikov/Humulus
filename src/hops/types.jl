# =============================================================================
# Hierarchy of Pure States (HOPS) type.
# =============================================================================

"""
    HOPS{N,MaxFockStates} <: Function

Functor implementing the hierarchy of pure states (HOPS).

The object stores reusable temporary buffers, a lookup table of square
roots used when applying creation and annihilation operators, and the
stochastic noise realization evaluated on a time grid.

# Type parameters
- `N`: number of effective modes.
- `MaxFockStates`: maximum number of retained Fock states per mode.

# Fields
- `f_tmp`: temporary storage for the values of the `f` functions.
- `g_tmp`: temporary storage for the values of the `g` functions.
- `_sqrt_of_int`: lookup table containing `sqrt.(0:MaxFockStates-1)`.
- `noise`: stochastic noise realization.
- `_time_grid`: time grid corresponding to `noise`.
"""
struct HOPS{N,MaxFockStates}<:Function
    f_tmp        :: MVector{N,ComplexF64}
    g_tmp        :: MVector{N,ComplexF64}
    _sqrt_of_int :: SVector{MaxFockStates,Float64}

    noise        :: Vector{ComplexF64}
    _time_grid   :: TimeGrid
end


"""
    HOPS{N,MaxFockStates}(time_grid)

Construct and return a `HOPS` functor for an `N`-mode system.

The `time_grid` determines the grid on which the stochastic noise
realization is stored.

`MaxFockStates` specifies the maximum number of retained Fock states per
mode.

# Arguments
- `time_grid::TimeGrid`: time grid used to store the stochastic noise realization.

# Exceptions

An exception is thrown if

- `N` or `MaxFockStates` is not an integer;
- `N ≤ 0`; or
- `MaxFockStates < 1`.
"""
function HOPS{N,MaxFockStates}(time_grid::TimeGrid) where {N,MaxFockStates}

    # Input validation. 
    N isa Int ||
        throw(ArgumentError("The number of modes `N` must be an `Int`, got $(typeof(N))."))

    MaxFockStates isa Int ||
        throw(ArgumentError("`MaxFockStates` must be an `Int`, got $(typeof(MaxFockStates))."))

    N > 0 ||
        throw(DomainError(N, "The number of modes `N` must be positive."))

    MaxFockStates >= 1 ||
        throw(DomainError(MaxFockStates, "`MaxFockStates` must be at least 1, since the vacuum state is always included."))

        
    # Allocate temporary buffers.
    f_tmp = MVector{N,ComplexF64}(undef)
    g_tmp = MVector{N,ComplexF64}(undef)

    # Precompute square roots of the occupation numbers.
    sqrt_of_int = SVector{MaxFockStates,Float64}(sqrt(i) for i in 0:MaxFockStates-1)

    grid_size = length(time_grid)
    noise     = zeros(ComplexF64, grid_size)

    # Construct the HOPS object.
    return HOPS{N,MaxFockStates}(f_tmp, g_tmp, sqrt_of_int, noise, time_grid)
end