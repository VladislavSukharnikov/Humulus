"""
    (sampler::NoiseSampler)(noise; checks=true)

Generate a complex Gaussian noise realization in place.

The sampler reuses an internal workspace to avoid repeated memory
allocations.

# Arguments
- `noise::Vector{ComplexF64}`: output vector. Its length must equal the
  discretization grid size.

# Keyword arguments
- `checks::Bool`: if `true`, verify that `noise` has the expected length.
"""
@inline function (sampler::NoiseSampler)(
                                    noise::Vector{ComplexF64}; 
                                    checks::Bool=true,
                                )

    (; _time_grid, _vecs, _vals, container) = sampler
    n_points = _time_grid.n_points

    checks && @assert length(noise) == n_points "Noise vector length mismatch."

    # Generate independent complex Gaussian coefficients.
    for i in 1:n_points
        container[i] = _vals[i] * (randn() + 1im * randn()) * invroot2
    end

    # Transform the independent Gaussian coefficients to correlated noise.
    mul!(noise, _vecs, container)

    return nothing
end