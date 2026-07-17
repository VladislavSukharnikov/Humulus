# -----------------------------------------------------------------------------
# Plotting mean Bloch vector components. 
# -----------------------------------------------------------------------------

using Plots

let ts_save=grid_params.ts_save, ρ_s=ρ_s
    σˣ, σʸ, σᶻ = bloch_vector(ρ_s)

    p1 = plot(ts_save, real(σˣ), ylims=(-1,1), xlabel="Time", label = "x");
    p2 = plot(ts_save, real(σʸ), ylims=(-1,1), xlabel="Time", label = "y");
    p3 = plot(ts_save, real(σᶻ), ylims=(-0.4,0), xlabel="Time", label = "z");

    plot(p1, p2, p3, layout=(3, 1), size=(700, 800), plot_title="Mean Bloch vector components.")
end