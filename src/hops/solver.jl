@inline function _solve_hops(
                        ::Val{N},
                        ::Val{MaxFockStates},
                        grid_params::GridParams,
                        solver_params::NamedTuple,
                        n_trajectories::Int,
                        path::String,
                    ) where {N, MaxFockStates}
    
    (; dt_max, ts_save) = grid_params

    # Create noise sampler.
    sampler!::NoiseSampler = sampler_from_cache(path, checks=false, clear_cache=false)::NoiseSampler;


    # Container for the physical density matrix
    ρ_s = zeros(ComplexF64, 2, 2, ts_save.n_points)
    traj_count = zeros(Int, 1)

        
    # Create HOPS object for the RHS evaluation.
    time_grid = sampler!._time_grid
    hops! = HOPS{N,MaxFockStates}(time_grid)

    
    # Create ODE problem.
    (; fock_dim, c_g, c_e) = solver_params
    u0    = init_hops(fock_dim, N, c_g, c_e)
    prob  = ODEProblem{true}(
                        ODEFunction{true}(hops!), 
                        u0, 
                        (ts_save.t_start, ts_save.t_end), 
                        solver_params,
                    )

    # Access noise vector from the problem's RHS functor. 
    noise = prob.f.f.noise

    for _ in 1:n_trajectories
        sampler!(noise)

        data = solve(
                prob, 
                Tsit5(), 
                adaptive=true, 
                dtmax=dt_max, 
                saveat=ts_save, 
                progress=true
            )

        traj_count[1] += 1

        for t_idx in 1:ts_save.n_points
            ψ = data.u[t_idx].x[1]
            save_data!(ρ_s, ψ, t_idx)
        end
        if traj_count[1]%10 == 0
            println("Finished $(traj_count[1]) trajectories.")
        end
    end

    return ArrayPartition(ρ_s, traj_count)
end


function solve_hops(
                grid_params::GridParams,
                bcf::BCF{N}, 
                atom_params::AtomParams,
                max_occupancies::Union{NTuple{N,Int},Int},
                n_trajectories::Int;
                noise_oversampling::Int = 2,
                KeyType::Type{<:Integer}=Int,
                IndType::Type{<:Integer}=Int,
                clear_cache::Bool=true,
            ) where {N}

    fock_space      = FockSpace(Val(N), max_occupancies, KeyType, IndType)
    solver_params   = create_solver_params(bcf, fock_space, atom_params)
    max_fock_states = maximum(max_occupancies)+1

    # Discretize the noise on a finer grid.
    (; dt_max, ts_save) = grid_params


    dt_noise  = dt_max / noise_oversampling
    grid_size = ceil(Int, (ts_save.t_end) / dt_noise)

    path = create_noise_cache(bcf, ts_save.t_end, grid_size)

    output = _solve_hops(
                    Val(N), 
                    Val(max_fock_states),  
                    grid_params, 
                    solver_params, 
                    n_trajectories,
                    path,
                )

    if clear_cache
        rm(path)
    end

    return output
end


function save_data!(ρ_s::Array{ComplexF64,3}, ψ::Array{ComplexF64,2}, t_idx::Int)
    ψ_g = ψ[1]
    ψ_e = ψ[2]
    
    C₀ = abs2(ψ_g)+abs2(ψ_e)

    ρ_s[1,1,t_idx] += abs2(ψ_g)/C₀
    ρ_s[1,2,t_idx] += conj(ψ_e)*ψ_g/C₀
    ρ_s[2,1,t_idx] += conj(ψ_g)*ψ_e/C₀
    ρ_s[2,2,t_idx] += abs2(ψ_e)/C₀
    
    return nothing
end
