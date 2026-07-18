"""
    TimeGrid <: AbstractVector{Float64}

Uniformly spaced one-dimensional time grid.

`TimeGrid` behaves like an `AbstractVector{Float64}` whose elements are
computed on demand rather than stored explicitly.

# Fields
- `t_start::Float64`: initial time.
- `t_end::Float64`: final time.
- `n_points::Int`: number of grid points.
"""
struct TimeGrid <: AbstractVector{Float64}
    t_start::Float64
    t_end::Float64
    n_points::Int

    function TimeGrid(t_start, t_end, n_points)
        n_points >= 2 ||
            throw(ArgumentError("`n_points` must be at least 2, got $n_points."))

        t_end >= t_start ||
            throw(ArgumentError("`t_end` must be greater than or equal to `t_start`."))

        new(Float64(t_start), Float64(t_end), Int(n_points))
    end
end

Base.size(grid::TimeGrid) = (grid.n_points,)

Base.IndexStyle(::Type{<:TimeGrid}) = IndexLinear()

function Base.getindex(grid::TimeGrid, i::Int)
    @boundscheck checkbounds(grid, i)

    # Uniform grid spacing.
    Δt = (grid.t_end - grid.t_start) / (grid.n_points - 1)

    # Compute the grid point on demand.
    return grid.t_start + (i - 1) * Δt
end