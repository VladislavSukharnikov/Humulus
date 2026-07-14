include("debug_utils.jl")

# =============================================================================
# Debugging: reconstruct the BCF from Monte Carlo noise sampling.
# =============================================================================
let 
    # Construct BCF.
    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)
 
    # Construct NoiseSampler.
    sampler! = let 
        t_end     = 20.0
        grid_size = 1000
        path = Humulus.create_noise_cache(bcf, t_end, grid_size)
        sampler! = Humulus.sampler_from_cache(
                                        path; 
                                        checks=true, 
                                        clear_cache=true,
                                    )
    end


    time_grid = sampler!._time_grid;
    grid_size = time_grid.n_points

    t_idx = rand(1:grid_size)

    # Number of realizations to sample.
    realization_num = 10_000;

    # Compute the reconstructed BCF with one time fixed, bcf(t, :), by averaging.
    reconstructed_bcf_t = average_noise!(sampler!, realization_num, t_idx);

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
end


# =============================================================================
# Debugging: Visually inspect the interpolation of the noise on the finer grid.
# =============================================================================
let 
    # Construct BCF.
    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)
 
    # Construct NoiseSampler.
    sampler! = let 
        t_end     = 20.0
        grid_size = 20
        path = Humulus.create_noise_cache(bcf, t_end, grid_size)
        sampler! = Humulus.sampler_from_cache(
                                        path; 
                                        checks=true, 
                                        clear_cache=true,
                                    )
    end

    # Extract the native sampling grid used by the noise sampler.
    coarse_time_grid = sampler!._time_grid
    coarse_grid_size = coarse_time_grid.n_points
    t_end            = coarse_time_grid[end]

    # Sample one noise trajectory on the native grid.
    coarse_noise = zeros(ComplexF64, coarse_grid_size)
    sampler!(coarse_noise)

    # Interpolate the sampled noise onto a finer grid.
    interpolation_grid_size = 1000
    interpolation_time_grid = range(0.0, t_end, interpolation_grid_size)

    interpolated_noise = zeros(ComplexF64, interpolation_grid_size)

    for t_idx in eachindex(interpolation_time_grid)
        interpolated_noise[t_idx] = Humulus.interpolate_noise(coarse_time_grid, coarse_noise, interpolation_time_grid[t_idx])
    end

    # Plot real part: coarse samples and interpolated curve.
    real_part_plot = scatter(coarse_time_grid, real.(coarse_noise), label = "Sampled noise")
    plot!(interpolation_time_grid, real.(interpolated_noise), label = "Interpolated noise")
    xlabel!("Time")
    ylabel!("Re[noise]")

    # Plot imaginary part: coarse samples and interpolated curve.
    imaginary_part_plot = scatter(coarse_time_grid, imag.(coarse_noise), label = "Sampled noise")
    plot!(interpolation_time_grid, imag.(interpolated_noise), label = "Interpolated noise")
    xlabel!("Time")
    ylabel!("Im[noise]")

    plot(real_part_plot, imaginary_part_plot, layout = (2, 1), size = (800, 700))
end


# =============================================================================
# Debugging: @code_warntype inspection of sampler!(z).
# =============================================================================
let 
    println("\n\n=============================================")
    println("Running @code_warntype for `sampler!(z)`")
    println("=============================================\n")

    # Construct BCF.
    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)
 
    # Construct NoiseSampler.
    sampler! = let 
        t_end     = 20.0
        grid_size = 1000
        path = Humulus.create_noise_cache(bcf, t_end, grid_size)
        sampler! = Humulus.sampler_from_cache(
                                        path; 
                                        checks=true, 
                                        clear_cache=true,
                                    )
    end
    
    grid_size = sampler!._time_grid.n_points
    z = zeros(ComplexF64, grid_size)

    # Trigger compilation. 
    sampler!(z)

    @code_warntype sampler!(z)
end


# =============================================================================
# Debugging: @code_warntype inspection of interpolate_noise(...).
# =============================================================================
let 
    println("\n\n===================================================")
    println("Running @code_warntype for `interpolate_noise(...)`")
    println("===================================================\n")

    # Construct BCF.
    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)
 
    # Construct NoiseSampler.
    sampler! = let 
        t_end     = 20.0
        grid_size = 1000
        path = Humulus.create_noise_cache(bcf, t_end, grid_size)
        sampler! = Humulus.sampler_from_cache(
                                        path; 
                                        checks=true, 
                                        clear_cache=true,
                                    )
    end
    
    time_grid = sampler!._time_grid
    grid_size = time_grid.n_points

    noise = zeros(ComplexF64, grid_size)
    sampler!(noise)

    # Initialize random time point.
    t = rand() * time_grid[end]

    # Trigger compilation. 
    Humulus.interpolate_noise(time_grid, noise, t)

    @code_warntype Humulus.interpolate_noise(time_grid, noise, t)
end


# =============================================================================
# Debugging: @code_warntype inspection of sampler_from_cache(...).
# =============================================================================
let 
    println("\n\n===================================================")
    println("Running @code_warntype for `sampler_from_cache(...)`")
    println("===================================================\n")

    # Construct BCF.
    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)
 
    # NoiseSampler construction options.
    t_end     = 20.0
    grid_size = 20

    # Create cache.
    path = Humulus.create_noise_cache(bcf, t_end, grid_size)

    # Trigger compilation.
    sampler! = Humulus.sampler_from_cache(
                        path, 
                        checks=true, 
                        clear_cache=true
                    )
    
    @code_warntype Humulus.sampler_from_cache(
                        path, 
                        checks=true, 
                        clear_cache=true
                    )
end


# =============================================================================
# Debugging: @code_warntype inspection of create_noise_cache(...).
# =============================================================================
let 
    println("\n\n===================================================")
    println("Running @code_warntype for `create_noise_cache(...)`")
    println("===================================================\n")

    # Construct BCF.
    N = 2
    bcf = Humulus.random_bcf(N)
 
    # Create noise cache.
    t_end   = 1.0;
    n_points = 100;

    # Trigger compilation.
    path = Humulus.create_noise_cache(bcf, t_end, n_points)
    rm(path)

    @code_warntype Humulus.create_noise_cache(bcf, t_end, n_points)
end