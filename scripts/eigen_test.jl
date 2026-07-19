using Humulus
using LinearAlgebra
using Plots
using BenchmarkTools

function average_noise!(
                    sampler!::Humulus.NoiseSampler{M},
                    trajectory_num::Int,
                    t_idx::Int
                ) where {M}

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

t_end = 20.0
grid_size = 2000
bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)

time_grid, covariance = Humulus._discretize_bcf(bcf, t_end, grid_size)

chol =cholesky!(covariance)
@time bcf_chol = Humulus.BCFCholesky(bcf, t_end, grid_size)

sampler1 = Humulus.NoiseSampler(bcf_chol.time_grid, bcf_chol.chol, zeros(ComplexF64, 2000))


time_grid = sampler1._time_grid;
    grid_size = time_grid.n_points

    t_idx = rand(1:grid_size)

    # Number of realizations to sample.
    realization_num = 10_000;

    # Compute the reconstructed BCF with one time fixed, bcf(t, :), by averaging.
    reconstructed_bcf_t = average_noise!(sampler1, realization_num, t_idx);

    # Compute the exact BCF, bcf(t, :).
    analytic_bcf_t = bcf.(time_grid[t_idx], time_grid)

    # Compare real parts of reconstructed and exact BCF.
    p_real = plot(time_grid, real.(analytic_bcf_t), label = "Exact BCF")
    plot!(time_grid, real.(reconstructed_bcf_t), label = "Reconstructed BCF")
    xlabel!("Time")
    ylabel!("Re[BCF]")

    # Compare imaginary parts of reconstructed and exact BCF.
    p_imag = plot(time_grid, imag.(analytic_bcf_t), label = "Exact BCF")
    plot!(time_grid, imag.(reconstructed_bcf_t), label = "Reconstructed BCF")
    xlabel!("Time")
    ylabel!("Im[BCF]")

    plot(p_real, p_imag, layout = (2, 1), size = (800, 700), plot_title="Exact vs reconstructed BCF.")