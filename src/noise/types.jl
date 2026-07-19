# =============================================================================
# Covariance Matrix Decompositions.
#
# The Cholesky decomposition is the default and recommended representation.
# The eigendecomposition is retained only for testing and validation.
# =============================================================================

abstract type AbstractBCFDecomposition end


"""
    BCFCholesky

Precomputed Cholesky factorization of the covariance matrix obtained by
discretizing the bath-correlation function (BCF).

The stored factor `chol` satisfies

    bcf_covariance = chol * chol'

and is used to generate correlated Gaussian noise efficiently.

# Fields
- `bcf`: bath-correlation function.
- `time_grid`: time grid used to discretize the BCF.
- `chol`: lower-triangular Cholesky factor of the covariance matrix.
"""
struct BCFCholesky <: AbstractBCFDecomposition
    # Bath-correlation function used to construct the matrix.
    bcf       :: BCF

    # Time grid on which the covariance matrix was discretized.
    time_grid :: TimeGrid

    # Cholesky decomposition.
    chol      :: LowerTriangular{ComplexF64, Matrix{ComplexF64}}
end


# =============================================================================
# Noise Generation. 
# =============================================================================

"""
    NoiseSampler

Generator of correlated Gaussian noise from a precomputed covariance-matrix
decomposition.

The sampler stores a matrix `A` satisfying

    bcf_covariance = A * A'

and generates correlated noise by multiplying a vector of independent complex
Gaussian random variables by `A`.

The matrix `A` may originate from any covariance-matrix decomposition
(e.g. Cholesky factorization or the eigendecomposition).

# Fields
- `_time_grid`: time grid used to discretize the covariance matrix.
- `_A`: matrix satisfying `bcf_covariance = A * A'`.
- `container`: reusable workspace for noise generation.
"""
struct NoiseSampler{M}
    # Time grid on which the covariance matrix was discretized.
    _time_grid::TimeGrid

    # Matrix A satisfying C = A*A'.
    _A::M

    # Storage reused during noise generation to avoid repeated allocations.
    container::Vector{ComplexF64}
end


# =============================================================================
# Legacy decomposition.
# =============================================================================

"""
    BCFEigen

Legacy covariance-matrix decomposition retained for testing and validation.

This type stores the scaled eigendecomposition of the covariance matrix,

    bcf_covariance = vecs * vecs'

where the columns of `vecs` are the eigenvectors scaled by the square roots
of the corresponding eigenvalues.

`BCFEigen` is **not** used by the production noise-generation pipeline.

# Fields
- `bcf`: bath-correlation function.
- `time_grid`: time grid used to discretize the BCF.
- `vecs`: eigenvectors scaled by the square roots of the corresponding
  eigenvalues.
"""
struct BCFEigen <: AbstractBCFDecomposition
    # Bath-correlation function used to construct the matrix.
    bcf       :: BCF

    # Time grid on which the covariance matrix was discretized.
    time_grid :: TimeGrid

    # Scaled eigenvectors of the covariance matrix.
    vecs      :: Matrix{ComplexF64}
end