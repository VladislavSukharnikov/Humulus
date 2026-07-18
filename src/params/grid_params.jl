# =============================================================================
# Grid parameters. 
# =============================================================================

"""
    GridParams

Time discretization parameters.

# Fields
- `ts_save::TimeGrid`: save-time grid.
- `dt_max::Float64`: maximum internal integrator step size.
"""
struct GridParams
    ts_save :: TimeGrid
    dt_max  :: Float64
end


"""
    GridParams(t_end, n_save, substeps)

Construct and return time discretization parameters.

The save-time grid is uniform on `[0, t_end]` with `n_save + 1` points.
The maximum internal step size is `Δt / substeps`, where `Δt` is the spacing
between consecutive save points.

# Arguments
- `t_end::Float64`: final simulation time.
- `n_save::Int`: number of save intervals.
- `substeps::Int`: maximum number of internal integrator substeps per save interval.

# Exceptions

An exception is thrown if

- `t_end ≤ 0`;
- `n_save < 1`; or
- `substeps < 1`.
"""
function GridParams(t_end::Float64, n_save::Int, substeps::Int)
    # Input validation.

    t_end > 0.0 ||
    throw(ArgumentError("`t_end` must be positive, got $t_end."))

    n_save >= 1 ||
        throw(ArgumentError("`n_save` must be at least 1, got $n_save."))

    substeps >= 1 ||
        throw(ArgumentError("`substeps` must be at least 1, got $substeps."))


    # Uniform save-time grid (includes t = 0).    
    ts_save = TimeGrid(0.0, t_end, n_save + 1)

    # Maximum internal integrator step size.
    dt_max = (ts_save[2] - ts_save[1])/substeps

    return GridParams(ts_save, dt_max)
end