"""
    BCFCholesky(bcf, t_end, grid_size)

Construct the Cholesky factorization of the covariance matrix obtained by
discretizing the bath correlation function (BCF).

The BCF is evaluated on a uniform time grid to assemble the corresponding
covariance matrix, which is then factorized using the Cholesky decomposition.
The resulting factorization can be reused to efficiently generate correlated
Gaussian noise realizations.

# Arguments
- `bcf::BCF`: bath correlation function.
- `t_end::Float64`: final time of the discretization interval.
- `grid_size::Int`: number of discretization points.

# Throws
- `ArgumentError`: if `grid_size < 1` or `t_end ≤ 0`.
- `PosDefException`: if the discretized covariance matrix is not positive
  definite.
"""
function BCFCholesky(bcf::BCF, t_end::Float64, grid_size::Int)
    # Construct the covariance matrix.
    time_grid, bcf_covariance = _discretize_bcf(bcf, t_end, grid_size)

    # Compute the Cholesky factorization of the covariance matrix.
    chol = cholesky!(bcf_covariance; check=true)

    # Construct and return BCFCholesky.
    return BCFCholesky(bcf, time_grid, chol.L)
end


"""
    _discretize_bcf(bcf, t_end, grid_size)

Discretize a bath-correlation function (BCF) on a uniform time grid.

The returned covariance matrix is stored as a `Hermitian` wrapper whose lower
triangular part contains the evaluated matrix entries.

# Arguments
- `bcf::BCF`: bath correlation function.
- `t_end::Float64`: final time of the discretization interval.
- `grid_size::Int`: number of discretization points.

# Returns
A tuple `(time_grid, bcf_covariance)`, where

- `time_grid::TimeGrid` is the discretization grid.
- `bcf_covariance::Hermitian{ComplexF64}` is the discretized covariance
  matrix.

# Throws
- `ArgumentError`: if `grid_size < 1`.
- `ArgumentError`: if `t_end ≤ 0`.
"""
function _discretize_bcf(bcf::BCF, t_end::Float64, grid_size::Int)    
    # Input validation.
    grid_size >= 2 ||
        throw(ArgumentError("`grid_size` must be at least 2, got $grid_size."))

    t_end > 0.0 ||
        throw(ArgumentError("`t_end` must be positive, got $t_end."))

    # Construct the discretization time grid.
    time_grid = TimeGrid(0.0, t_end, grid_size)

    # Assemble the lower triangular part of the covariance matrix.
    bcf_covariance = Matrix{ComplexF64}(undef, grid_size, grid_size)
    @inbounds for t_idx in 1:grid_size
        t = time_grid[t_idx]
        for s_idx in 1:t_idx
            bcf_covariance[t_idx, s_idx] = bcf(t, time_grid[s_idx])
        end
    end

    return time_grid, Hermitian(bcf_covariance, :L)
end



# =============================================================================
# Legacy eigendecomposition implementation.
#
# Retained for testing and validation. The production implementation uses the
# Cholesky factorization exclusively.
# =============================================================================

"""
    BCFEigen(bcf, t_end, grid_size)

Construct the legacy eigendecomposition of the covariance matrix obtained by
discretizing the bath-correlation function (BCF).

This implementation is retained for testing and validation. The production
noise-generation pipeline uses `BCFCholesky`, which is substantially more
efficient while producing equivalent correlated Gaussian noise realizations.

# Arguments
- `bcf::BCF`: bath-correlation function.
- `t_end::Float64`: final time of the discretization interval.
- `grid_size::Int`: number of discretization points.
"""
function BCFEigen(bcf::BCF, t_end::Float64, grid_size::Int)
    # Construct the covariance matrix.
    time_grid, bcf_covariance = _discretize_bcf(bcf, t_end, grid_size)

    # Compute the eigendecomposition of the covariance matrix.
    vals, vecs = eigen!(bcf_covariance)

    # Warn if the matrix is not positive semidefinite.
    any(<(0.0), vals) && @warn "Negative eigenvalue detected."

    # Scale each eigenvector by the square root of its corresponding eigenvalue.
    vecs .*= reshape(sqrt.(real(vals)), 1, :)

    # Construct and return BCFEigen.
    return BCFEigen(bcf, time_grid, vecs)
end