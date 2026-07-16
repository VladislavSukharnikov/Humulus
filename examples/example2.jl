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
    clear_cache=false,
    show_progress=true,
);

# -----------------------------------------------------------------------------
# Sample trajectories.
# -----------------------------------------------------------------------------

# # # # # # # # # # #

using Distributed

addprocs(5)

@everywhere using Humulus

n_trajectories = 100
out = @time solve_hops(
    grid_params,
    bcf,
    atom_params,
    max_occupancy,
    n_trajectories;
    clear_cache=false,
    show_progress=true,
    workers=workers(),
);

ρ_s = out.x[1]./out.x[2]

@info "Expected number of trajectories: $(n_trajectories*length(workers())) "
@info "Total number of trajectories: $(out.x[2][1])"