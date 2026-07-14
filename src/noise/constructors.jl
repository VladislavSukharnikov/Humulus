"""
    BCFEigen(bcf, t_end, grid_size)

Construct a `BCFEigen` object from a bath-correlation function (BCF).

The BCF is discretized on a uniform time grid, and the eigendecomposition
of the resulting kernel matrix is computed. The returned object stores the
time grid together with the eigendecomposition, which can be saved and
later reused for efficient noise generation.

# Arguments
- `bcf`: bath-correlation function.
- `t_end`: final time of the discretization interval.
- `grid_size`: number of discretization points.

# Returns
- `BCFEigen`: precomputed eigendecomposition of the discretized BCF.
"""
function BCFEigen(
                bcf::BCF, 
                t_end::Float64, 
                grid_size::Int,
            )
        
    # Construct the discretization grid.
    time_grid = TimeGrid(0.0, t_end, grid_size)

    # Assemble the discretized BCF kernel matrix.
    bcf_matrix = zeros(ComplexF64, grid_size, grid_size)
    for t in 1:grid_size, s in 1:grid_size
        bcf_matrix[t,s] = bcf(time_grid[t], time_grid[s])
    end

    # Compute the eigendecomposition of the BCF kernel.
    vals, vecs = eigen(bcf_matrix)

    # Warn if the kernel is not positive semidefinite.
    if any(real(vals) .< 0)
        @warn "Negative eigenvalue detected."
    end

    return BCFEigen(bcf, time_grid, sqrt.(real(vals)), vecs)
end