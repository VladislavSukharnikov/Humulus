struct HOPS{N,MaxFockStates}<:Function
    f_tmp        :: MVector{N,ComplexF64}
    g_tmp        :: MVector{N,ComplexF64}
    _sqrt_of_int :: SVector{MaxFockStates,Float64}

    noise        :: Vector{ComplexF64}
    _time_grid   :: TimeGrid
end

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
    f_tmp = zeros(MVector{N,ComplexF64})
    g_tmp = zeros(MVector{N,ComplexF64})

    # Precompute square roots of the occupation numbers.
    sqrt_of_int  = SVector{MaxFockStates,Float64}(sqrt.(0:MaxFockStates-1))

    grid_size = length(time_grid)
    noise     = zeros(ComplexF64, grid_size)

    # Construct the HOPS object.
    return HOPS{N,MaxFockStates}(f_tmp, g_tmp, sqrt_of_int, noise, time_grid)
end
