"""
    BCFEigen(bcf, t_end, grid_size)

Construct and return the eigendecomposition of a discretized bath
correlation function (BCF).

The BCF is evaluated on a uniform time grid to form the corresponding bath
correlation matrix. The returned object stores the time grid together with
the eigendecomposition of this matrix, which can be saved and later reused
for efficient noise generation.

# Arguments
- `bcf::BCF`: bath correlation function.
- `t_end::Float64`: final time of the discretization interval.
- `grid_size::Int`: number of discretization points.
"""
function BCFEigen(
                bcf::BCF, 
                t_end::Float64, 
                grid_size::Int,
            )

    # Input validation.    
    grid_size >= 1 ||
        throw(ArgumentError("`grid_size` must be at least 1, got $grid_size."))

    t_end > 0.0 ||
        throw(ArgumentError("`t_end` must be positive, got $t_end."))

    # Construct the discretization grid.
    time_grid = TimeGrid(0.0, t_end, grid_size)

    # Assemble the upper triangular part of the discretized bath correlation matrix.
    bcf_matrix = Matrix{ComplexF64}(undef, grid_size, grid_size)
    @inbounds for t in 1:grid_size
        for s in 1:t
            bcf_matrix[t, s] = bcf(time_grid[t], time_grid[s])
        end
    end

    vals, vecs = eigen!(Hermitian(bcf_matrix, :L))

    # Warn if the matrix is not positive semidefinite.
    if any(vals .< 0.0)
        @warn "Negative eigenvalue detected."
    end

    return BCFEigen(bcf, time_grid, sqrt.(vals), vecs)
end