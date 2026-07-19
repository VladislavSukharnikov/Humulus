# This example reproduces results from Fig. 1 of
# Physical Review Research 8, 013123 (2026)
# using the hierarchy of master equations (HME).

using Humulus

# -----------------------------------------------------------------------------
# Squeezed reservoir parameters
# -----------------------------------------------------------------------------

bcf = let
    ω₀ = 5.0       # central (carrier) frequency
    Γ  = 1.0       # spectral half-width
    r  = 1.5       # squeezing parameter
    φ  = 0.0       # squeezing phase
    γ  = 1.0       # atom–field coupling strength
    
    one_mode_squeezed_bcf(ω₀, Γ, r, φ, γ)
end;

# bcf = let
#     ω₀ = 5.0        # central (carrier) frequency
#     Γ₀ = 2.0        # input spectral width
#     κ  = 1.0        # cavity loss rate
#     ϵ  = 0.5        # effective DPA pump amplitude
#     φ  = Float64(π) # squeezing phase
#     γ  = 1.0        # atom–field coupling strength

#     bcf = three_mode_squeezed_bcf(ω₀, Γ₀, κ, ϵ, φ, γ)
# end;

# -----------------------------------------------------------------------------
# Atomic parameters
# -----------------------------------------------------------------------------

atom_params = let
    # central (carrier) frequency
    ν₀ = 5.0

    # Initial atomic state:
    # |ψ⟩ = c_g |g⟩ + c_e |e⟩
    c_e = inv(√2)
    c_g = inv(√2) * exp(-1im * π/4)

    AtomParams(ν₀, c_g, c_e)
end;


# -----------------------------------------------------------------------------
# Time discretization and propagation
# -----------------------------------------------------------------------------

grid_params = let
    t_end    = 20.0  # final simulation time
    n_save   = 250   # number of intervals between saved time points
    substeps = 5     # maximum internal integration substeps per save interval
    @info "Expected minimal number of steps is $(n_save*substeps)."
    GridParams(t_end, n_save, substeps)
end;


# -----------------------------------------------------------------------------
# Pseudo-Fock-space truncation
# -----------------------------------------------------------------------------

# truncation of pseudo-Fock basis. 
max_occupancy = 50


# -----------------------------------------------------------------------------
# Propagate the hierarchy of master equations.
# -----------------------------------------------------------------------------

ρ_s = @time solve_hme(
    grid_params,
    bcf,
    atom_params,
    max_occupancy,
);


# -----------------------------------------------------------------------------
# Running HOPS on the main process. 
# -----------------------------------------------------------------------------

n_trajectories = 1
out = @time solve_hops(
    grid_params,
    bcf,
    atom_params,
    max_occupancy;
    n_trajectories=n_trajectories,
    clear_cache=false,
    show_progress=true,
    logging=true,
);

ρ_s = out.x[1]./out.x[2];

# -----------------------------------------------------------------------------
# Distributed computation of HOPS.
# -----------------------------------------------------------------------------

using Distributed

n_workers = Sys.CPU_THREADS
addprocs(n_workers)

@everywhere using Humulus

n_trajectories = 100
n_batches = 100
@info "Expected number of trajectories: $(n_trajectories*length(workers())*n_batches) "

out = @batched n_batches solve_hops(
                            grid_params,
                            bcf,
                            atom_params,
                            max_occupancy;
                            n_trajectories=n_trajectories,
                            clear_cache=false,
                            show_progress=false,
                            logging=false,
                            workers=workers(),
                        );



ρ_s = out.x[1]./out.x[2]

@info "Actual number of trajectories: $(out.x[2][1])"