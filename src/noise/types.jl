# =============================================================================
# Noise Parameters.
# =============================================================================

"""
    BCFEigen

Precomputed eigendecomposition of a discretized bath-correlation function (BCF).

# Fields
- `bcf`: bath-correlation function.
- `time_grid`: time grid used to discretize the BCF.
- `vals`: square roots of the eigenvalues of the discretized BCF kernel.
- `vecs`: eigenvectors of the discretized BCF kernel.
"""
struct BCFEigen
    # Bath-correlation function used to construct the matrix.
    bcf::BCF

    # Time grid on which the BCF was discretized.
    time_grid::TimeGrid

    # Square roots of the eigenvalues of the discretized BCF kernel.
    vals::Vector{Float64}

    # Eigenvectors of the discretized BCF kernel.
    vecs::Matrix{ComplexF64}
end


# =============================================================================
# Noise Generation. 
# =============================================================================
"""
    NoiseSampler

Generator of stochastic noise realizations from a precomputed bath
correlation function (BCF) eigendecomposition.

The type stores the discretization grid, the eigendecomposition of the
discretized BCF kernel, and a reusable workspace to avoid repeated memory
allocations during noise generation.

# Fields
- `_time_grid`: time grid used to discretize the BCF.
- `_vals`: square roots of the eigenvalues of the discretized BCF kernel.
- `_vecs`: eigenvectors of the discretized BCF kernel.
- `container`: reusable workspace for noise generation.
"""
struct NoiseSampler
    # Time grid on which the BCF was discretized.
    _time_grid::TimeGrid

    # Square roots of the eigenvalues of the discretized BCF kernel.
    _vals::Vector{Float64}

    # Eigenvectors of the discretized BCF kernel.
    _vecs::Matrix{ComplexF64}

    # Storage reused during noise generation to avoid repeated allocations.
    container::Vector{ComplexF64}
end