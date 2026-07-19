"""
    (sampler::NoiseSampler)(noise; checks=true)

Generate one realization of correlated complex Gaussian noise in-place.

The sampler first generates independent standard complex Gaussian random
variables and then applies the stored covariance-matrix factor to obtain the
desired correlations. The internal workspace is reused to avoid repeated
memory allocations.

# Arguments
- `noise::Vector{ComplexF64}`: output vector. Its length must equal the
  discretization grid size.

# Keyword arguments
- `checks::Bool=true`: if `true`, verify that `noise` has the expected length.
"""
@inline function (sampler::NoiseSampler)(
                                    noise::Vector{ComplexF64}, 
                                    checks::Bool=true,
                                )

    (; _time_grid, _A, container) = sampler
    n_points = _time_grid.n_points

    # Parameter validation.
    checks && @assert length(noise)==length(container)==n_points "Noise vector length mismatch."
    checks && @assert size(_A) == (n_points, n_points) "Expected a $n_points × $n_points covariance factor, got $(size(_A))."

    # Generate independent complex Gaussian coefficients.
    @inbounds for i in 1:n_points
        container[i] = (randn() + 1im * randn()) * invroot2
    end

    # Transform the independent Gaussian variables to correlated noise.
    mul!(noise, _A, container)

    return nothing
end