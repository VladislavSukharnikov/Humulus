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
end


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
end


# -----------------------------------------------------------------------------
# Time discretization and propagation
# -----------------------------------------------------------------------------

grid_params = let
    t_end    = 20.0  # final simulation time
    n_save   = 250   # number of intervals between saved time points
    substeps = 4     # maximum internal integration substeps per save interval
    @info "Expected minimal number of steps is $(n_save*substeps)."
    GridParams(t_end, n_save, substeps)
end


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
# Plotting mean Bloch vector components. 
# -----------------------------------------------------------------------------

using Plots

let ts_save=grid_params.ts_save
    σˣ, σʸ, σᶻ = bloch_vector(ρ_s)

    p1 = plot(ts_save, real(σˣ), ylims=(-1,1), xlabel="Time", label = "x");
    p2 = plot(ts_save, real(σʸ), ylims=(-1,1), xlabel="Time", label = "y");
    p3 = plot(ts_save, real(σᶻ), ylims=(-0.4,0), xlabel="Time", label = "z");

    plot(p1, p2, p3, layout=(3, 1), size=(700, 800), plot_title="Mean Bloch vector components.")
end



# # # # # # # # # # #



n_trajectories = 1000
out = @time solve_hops(
    grid_params,
    bcf,
    atom_params,
    max_occupancy,
    n_trajectories;
    clear_cache=false,
);

ρ_s_av = out.x[1] ./ out.x[2]

let ts_save = grid_params.ts_save;
    σˣ, σʸ, σᶻ = bloch_vector(ρ_s_av);

    p1 = plot(ts_save, real(σˣ), ylims=(-1,1), xlabel="Time", title="Mean x-component", label = "x");
    p2 = plot(ts_save, real(σʸ), ylims=(-1,1), xlabel="Time", title="Mean y-component", label = "y");
    p3 = plot(ts_save, real(σᶻ), ylims=(-0.4,0), xlabel="Time", title="Mean z-component", label = "z");

    plot(p1, p2, p3, layout=(3, 1), size=(700, 800))
end


max_occupancies = 50
fock_space      = Humulus.FockSpace(Val(N), max_occupancies, Int, Int)
    solver_params   = Humulus.create_solver_params(bcf, fock_space, atom_params)
    max_fock_states = maximum(max_occupancies)+1

(; ts_save, dt_max) = grid_params;

    # Container for the physical density matrix
    ρ_s = zeros(ComplexF64, 2, 2, ts_save.n_points)

    (; fock_dim, c_g, c_e) = solver_params
    hme! = Humulus.HME{N,MaxFockStates}();
    ρ_0  = Humulus.init_hme(fock_dim, c_g, c_e)
    prob = ODEProblem{true}(ODEFunction{true}(hme!), ρ_0, (ts_save.t_start, ts_save.t_end), solver_params);

    data = solve(
            prob, 
            Tsit5(), 
            adaptive=true, 
            dtmax=dt_max, 
            saveat=ts_save, 
            dense = false, 
            save_everystep = false, 
            calck = false,
            progress = true,
        )