
"""
    solve_hops(grid_params, bcf, atom_params, max_occupancies;
               n_trajectories=1,
               noise_oversampling=2,
               KeyType=Int,
               IndType=Int,
               clear_cache=true,
               show_progress=true,
               workers=[1])

Solve the hierarchy of pure states (HOPS).

This method constructs the pseudo-Fock space and solver parameters before
integrating the HOPS trajectories.

# Arguments
- `grid_params::GridParams`: time discretization parameters.
- `bcf::BCF`: bath correlation function.
- `atom_params::AtomParams`: atomic parameters and initial state.
- `max_occupancies::Union{NTuple{N,Int},Int,Vector{Int}}`: pseudo-Fock-space truncation.

# Keyword arguments
- `n_trajectories::Int`: number of stochastic trajectories.
- `noise_oversampling::Int`: oversampling factor used to discretize the noise.
- `KeyType::Type{<:Integer}`: integer type used for pseudo-Fock-space keys.
- `IndType::Type{<:Integer}`: integer type used for pseudo-Fock-space indices.
- `clear_cache::Bool`: remove the cached BCF eigenvalue decomposition after the computation.
- `show_progress::Bool`: display a progress bar during the computation.
- `workers::Vector{Int}`: worker processes used for trajectory parallelization.

# Exceptions

Any exceptions thrown while constructing the pseudo-Fock space, generating
the stochastic noise, or integrating the HOPS trajectories are propagated.
"""
function solve_hops(
                grid_params::GridParams,
                bcf::BCF{N}, 
                atom_params::AtomParams,
                max_occupancies::Union{NTuple{N,Int},Int};
                n_trajectories::Int=1,
                noise_oversampling::Int=4,
                KeyType::Type{<:Integer}=Int,
                IndType::Type{<:Integer}=Int,
                clear_cache::Bool=true,
                show_progress::Bool=true,
                logging::Bool=true,
                workers::Vector{Int}=[1],
            ) where {N}

    # Input validation.    
    n_trajectories ≥ 0 || 
        throw(ArgumentError("n_trajectories must be non-negative, got $n_trajectories"))

    noise_oversampling ≥ 2 || 
        throw(ArgumentError("noise_oversampling must be at least 2, got $noise_oversampling"))

    # Construct the pseudo-Fock space and solver parameters.
    fock_space      = FockSpace(Val(N), max_occupancies, KeyType, IndType)
    solver_params   = create_solver_params(bcf, fock_space, atom_params)
    max_fock_states = maximum(max_occupancies)+1

    # Construct the grid used to sample the stochastic noise.
    (; dt_max, ts_save) = grid_params
    dt_noise  = dt_max / noise_oversampling
    grid_size = ceil(Int, (ts_save.t_end) / dt_noise)

    # Compute or retrieve the cached BCF eigenvalue decomposition.
    path = get_cache(BCFCholesky, bcf, ts_save.t_end, grid_size; logging=logging)

    # Solve the HOPS equations serially or in parallel.
    if workers==[1] || length(workers)==1 || length(workers)==0
        output = _solve_hops(
                        Val(N), 
                        Val(max_fock_states),  
                        grid_params, 
                        solver_params, 
                        n_trajectories,
                        path,
                        show_progress,
                    )
    else
        wp = WorkerPool(workers)
        jobs = fill(n_trajectories, length(workers))
        data_all = pmap(wp, jobs) do ntraj
                        _solve_hops(
                            Val(N),
                            Val(max_fock_states),
                            grid_params,
                            solver_params,
                            ntraj,
                            path,
                            show_progress,
                            )
                end
        output = sum(data_all[:])
    end

    # Remove the cached decomposition if requested.
    clear_cache && rm(path)

    return output
end

# Internal implementation specialized on the number of modes and the maximum
# number of retained Fock states.
@inline function _solve_hops(
                        ::Val{N},
                        ::Val{MaxFockStates},
                        grid_params::GridParams,
                        solver_params::NamedTuple,
                        n_trajectories::Int,
                        path::String,
                        show_progress::Bool,
                    ) where {N,MaxFockStates}
    (; dt_max, ts_save) = grid_params

    # Construct the noise sampler.
    sampler!::NoiseSampler = sampler_from_cache(path; checks=false, clear_cache=false)::NoiseSampler

    # Allocate the output density matrix.
    ρ_s = zeros(ComplexF64, 2, 2, ts_save.n_points)
        
    # Create HOPS object for the RHS evaluation.
    time_grid = sampler!._time_grid
    hops! = HOPS{N,MaxFockStates}(time_grid)

    # Create ODE problem.
    (; fock_dim, c_g, c_e) = solver_params
    u0 = init_hops(fock_dim, N, c_g, c_e)
    prob = ODEProblem{true}(
                        ODEFunction{true}(hops!), 
                        u0, 
                        (ts_save.t_start, ts_save.t_end), 
                        solver_params,
                    )

    # Access the noise vector stored by the HOPS functor.
    noise = prob.f.f.noise

    # Accumulate stochastic trajectories.
    traj_count = 0
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

        traj_count += 1

        for t_idx in 1:ts_save.n_points
            ψ = data.u[t_idx].x[1]
            accumulate_trajectory!(ρ_s, ψ, t_idx)
        end
        
        # Displaying the progress.
        if show_progress
            step = max(1, n_trajectories ÷ 10)  # message every ~10%
            if traj_count % step == 0 || traj_count == n_trajectories
                percent = round(Int, 100 * traj_count / n_trajectories)
                @info "Finished $percent% of trajectories."
            end
        end
    end

    return ArrayPartition(ρ_s, [traj_count])
end


"""
    accumulate_trajectory!(ρ_s, ψ, t_idx)

Accumulate the normalized reduced density matrix corresponding to the
stochastic state `ψ` into `ρ_s[:, :, t_idx]`.
"""
function accumulate_trajectory!(ρ_s::Array{ComplexF64,3}, ψ::Array{ComplexF64,2}, t_idx::Int)
    ψ_g = ψ[1]
    ψ_e = ψ[2]

    # Normalize the stochastic state before accumulating its contribution.
    C₀ = abs2(ψ_g)+abs2(ψ_e)

    ρ_s[1,1,t_idx] += abs2(ψ_g)/C₀
    ρ_s[1,2,t_idx] += conj(ψ_e)*ψ_g/C₀
    ρ_s[2,1,t_idx] += conj(ψ_g)*ψ_e/C₀
    ρ_s[2,2,t_idx] += abs2(ψ_e)/C₀
    
    return nothing
end


"""
    @batched n_batches expr

Evaluate `expr` repeatedly for `n_batches` iterations, accumulating the
returned values.

# Arguments
- `n_batches`: number of batches.
- `expr`: expression to evaluate in each batch.
"""
macro batched(n_batches, call)
    quote
        local _n_batches = $(esc(n_batches))

        # Accumulated result from all completed batches.
        local _result = nothing

        # Number of successfully completed batches.
        local _completed = 0

        try
            for _batch in 1:_n_batches
                # Run the next batch.
                local _out = $(esc(call))

                # Accumulate the batch result.
                if _result === nothing
                    _result = _out
                else
                    _result += _out
                end

                _completed = _batch
                @info "Finished batch $_batch/$_n_batches."
            end
        catch err
            @warn "Stopped after $_completed completed batches." exception=(err, catch_backtrace())
        end

        _result
    end
end