# This example reproduces plots from Fig. 4 (b-c) of
# PHYSICAL REVIEW RESEARCH 8, 013123 (2026)
# by using HOPS.

using Humulus


# -----------------------------------------------------------------------------
# Squeezed reservoir parameters
# -----------------------------------------------------------------------------

bcf = let
    ω₀ = 5.0        # central (carrier) frequency
    Γ₀ = 2.0        # input spectral width
    κ  = 1.0        # cavity loss rate
    ϵ  = 0.5        # effective DPA pump amplitude
    φ  = Float64(π) # squeezing phase
    γ  = 1.0        # atom–field coupling strength

    bcf = three_mode_squeezed_bcf(ω₀, Γ₀, κ, ϵ, φ, γ)
end


# -----------------------------------------------------------------------------
# Atomic parameters
# -----------------------------------------------------------------------------

atom_params = let
    # central (carrier) frequency
    ω₀ = 5.0

    # Initial atomic state:
    # |ψ⟩ = c_g |g⟩ + c_e |e⟩
    c_e = inv(√2)
    c_g = inv(√2) * exp(-1im * π / 4)

    AtomParams(ω₀, c_g, c_e)
end


# -----------------------------------------------------------------------------
# Time discretization and propagation
# -----------------------------------------------------------------------------

grid_params = let
    t_end    = 20.0   # final simulation time
    n_save   = 500    # number of intervals between saved time points
    substeps = 2     # maximum internal integration substeps per save interval
    @info "Expected minimal number of steps is $(n_save*substeps)."
    GridParams(t_end, n_save, substeps)
end


# -----------------------------------------------------------------------------
# Pseudo-Fock-space truncation
# -----------------------------------------------------------------------------

# maximum sum of all occupation numbers. 
max_occupancy = 10


# -----------------------------------------------------------------------------
# Sample trajectories.
# -----------------------------------------------------------------------------

n_trajectories = 1
out = @time solve_hops(
    grid_params,
    bcf,
    atom_params,
    max_occupancy,
    n_trajectories,
    clear_cache=false
);

# -----------------------------------------------------------------------------
# Sample trajectories.
# -----------------------------------------------------------------------------

using Distributed


ρ_s = out.x[1] ./ out.x[2]

job_list = [1000 for _ in 1:workers_num]

wp = WorkerPool(workers())  # workers is Vector{Int} of worker IDs


out = pmap(wp, job_list) do n_traj
        integrate_trajectories(path, n_traj, bcf, grid_params, atom_params, max_occupancies)
    end

    out_sum = sum(out)

    
# -----------------------------------------------------------------------------
# Plotting mean Bloch vector components
# -----------------------------------------------------------------------------

using Plots

let ts_save=grid_params.ts_save, ρ_s=ρ_s
    σˣ, σʸ, σᶻ = bloch_vector(ρ_s, ts_save.n_points)

    p1 = plot(ts_save, real(σˣ), ylims=(-1,1), xlabel="Time", title="Mean x-component", label = "x");
    p2 = plot(ts_save, real(σʸ), ylims=(-1,1), xlabel="Time", title="Mean y-component", label = "y");
    p3 = plot(ts_save, real(σᶻ), ylims=(-0.4,0), xlabel="Time", title="Mean z-component", label = "z");

    plot(p1, p2, p3, layout=(3, 1), size=(700, 800))
end