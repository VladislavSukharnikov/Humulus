struct TimeGrid <: AbstractVector{Float64}
    t_start::Float64
    t_end::Float64
    n_points::Int

    function TimeGrid(t_start, t_end, n_points)
        n_points >= 2 || throw(ArgumentError("n_points must be at least 2"))
        t_end >= t_start || throw(ArgumentError("t_end must be >= t_start"))

        new(Float64(t_start), Float64(t_end), Int(n_points))
    end
end

Base.size(grid::TimeGrid) = (grid.n_points,)
Base.IndexStyle(::Type{<:TimeGrid}) = IndexLinear()

function Base.getindex(grid::TimeGrid, i::Int)
    @boundscheck checkbounds(grid, i)

    Δt = (grid.t_end - grid.t_start) / (grid.n_points - 1)
    return grid.t_start + (i - 1) * Δt
end