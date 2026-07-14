struct HOPS{N,MaxFockStates}<:Function
    f_tmp        :: MVector{N,ComplexF64}
    g_tmp        :: MVector{N,ComplexF64}
    _sqrt_of_int :: SVector{MaxFockStates,Float64}

    noise        :: Vector{ComplexF64}
    _time_grid   :: TimeGrid
end

function HOPS{N,MaxFockStates}(time_grid::TimeGrid) where {N,MaxFockStates}
    @assert N isa Int "Number of modes must be integer."
    @assert N > 0     "Number of modes must be positive."

    @assert MaxFockStates isa Int "`MaxFockStates` must be integer."
    @assert MaxFockStates > 1     "`MaxFockStates` must be larger or equal 1."

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
