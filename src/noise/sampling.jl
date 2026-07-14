"""
    (sample_noise::NoiseSampler)(noise)

Generate a complex Gaussian noise realization in-place.

Internal workspace is reused to avoid repeated memory allocations.

# Arguments
- `noise`: output vector. Its length must equal the discretization grid size.
"""
@inline function (sampler::NoiseSampler)(
                                    noise::Vector{ComplexF64}; 
                                    checks::Bool=true,
                                )

    (; _time_grid, _vecs, _vals, container) = sampler

    n_points = _time_grid.n_points

    if checks
        @assert length(noise) == n_points "Noise vector length mismatch."
    end

    # Generate independent complex Gaussian coefficients.
    for i in 1:n_points
        container[i] = _vals[i] * (randn() + 1im * randn()) * invroot2
    end

    mul!(noise, _vecs, container)

    return nothing
end