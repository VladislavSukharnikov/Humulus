using Humulus
using JLD2
using RecursiveArrayTools
using Plots

# =============================================================================
# Utilities.
# =============================================================================

function kwwrapper(mod::Module, f::Function)
    name = string(nameof(f))
    sym = only(filter(s ->
        occursin("#$name#", String(s)) &&
        !startswith(String(s), "##"),
        names(mod; all=true)))
    getfield(mod, sym)
end

"""
    average_noise!(sampler!, trajectory_num, t_idx)

Estimate the bath-correlation function by Monte Carlo averaging.

Generates `trajectory_num` independent noise realizations and returns the
sample estimate of

z[t_idx] * conj(z[:]⟩

The returned vector therefore represents the noise autocorrelation 
function between the fixed time index `t_idx` and all grid times.

# Arguments
- `sampler!::NoiseSampler`: Callable noise sampler.
- `trajectory_num::Int`: Number of independent noise realizations.
- `t_idx::Int`: Index of the fixed time point.

# Returns
- `Vector{ComplexF64}`: Monte Carlo estimate of the BCF evaluated between
  `t_idx` and every grid point.
"""
function average_noise!(
                    sampler!::Humulus.NoiseSampler,
                    trajectory_num::Int,
                    t_idx::Int
                )

    grid_size = sampler!._time_grid.n_points

    # Storage for a single noise realization.
    z = zeros(ComplexF64, grid_size)

    # Accumulator for the Monte Carlo estimate
    bcf_t = zeros(ComplexF64, grid_size)

    for _ in 1:trajectory_num
        # Generate a single noise realization.
        sampler!(z)

        # Accumulate the correlation between the fixed time point `t_idx`
        # and every grid point.
        @. bcf_t += z[t_idx] * conj(z)
    end

    # Compute the Monte Carlo average.
    return bcf_t ./ trajectory_num
end