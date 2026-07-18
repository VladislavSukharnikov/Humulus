"""
    (hops::HOPS{N})(du, u, solver_params, t)

Evaluate the right-hand side of the hierarchy of pure states (HOPS).

The result is written in-place to `du`, overwriting its previous contents.

# Arguments
- `du`: output state derivative.
- `u`: hierarchy state.
- `solver_params`: solver parameters created by `create_solver_params`.
- `t::Float64`: time.
"""
function (hops::HOPS{N})(du::ArrayPartition{ComplexF64},
                           u::ArrayPartition{ComplexF64}, 
                           solver_params::NamedTuple, 
                           t::Float64) where {N}
    (; Γ, G, f_vector, g_vector) = solver_params
    (; f_tmp, g_tmp, noise) = hops

    ψ, X    = u.x
    dψ, dX  = du.x

    # Evaluate the average dipole moment.
    L_av :: ComplexF64 = (conj(ψ[1])*ψ[2]+conj(ψ[2])*ψ[1])/(abs2(ψ[1])+abs2(ψ[2]))

    # Evaluate the noise value.
    z_t = interpolate_noise(hops._time_grid, noise, t)

    # Evaluate the auxiliary quantities.
    for j in 1:N
        f_tmp[j] = 1im*G[j]*f_vector[j](t)
        g_tmp[j] = 1im*G[j]*conj(g_vector[j](t))
        z_t     -= 1im*f_vector[j](t)*X[j]
        dX[j]    = -Γ[j]*X[j]+G[j]^2*conj(g_vector[j](t))*L_av
    end

    # Compute the HOPS right-hand side.
    _fill_dψ!(dψ, ψ, solver_params, L_av, z_t, f_tmp, g_tmp, hops._sqrt_of_int) 

    return nothing
end


"""
    _fill_dψ!(dψ, ψ, solver_params, L_av, z_t, f_tmp, g_tmp, sqrt_of_int)

Compute the HOPS state derivative in place.

The result is written to `dψ`, overwriting its previous contents.
"""
@inline function _fill_dψ!(dψ::AbstractArray{ComplexF64,2},
                           ψ::AbstractArray{ComplexF64,2},
                           solver_params::NamedTuple, 
                           L_av::ComplexF64,
                           z_t::ComplexF64,
                           f_tmp::MVector,
                           g_tmp::MVector,
                           sqrt_of_int::SVector)
    (; N, ν₀, basis_states, fock_dim, Γ, raise_index, lower_index) = solver_params
    mode_range = 1:N

    @inbounds for n_idx in 1:fock_dim
        n = basis_states[n_idx]
        ψ₁, ψ₂ = @view ψ[:,n_idx]

        dψ₁ = 1im*(ν₀/2*ψ₁-conj(z_t)*ψ₂)
        dψ₂ =-1im*(ν₀/2*ψ₂+conj(z_t)*ψ₁)

        totalΓ = 0.0
    
        # sqrt_of_int[k] = √(k-1).
        # Accumulate couplings to neighboring hierarchy elements
        @inbounds for mode in mode_range
            nₐ = n[mode]; 
            totalΓ += Γ[mode]*nₐ
            idx_prev, idx_next = lower_index[mode,n_idx], raise_index[mode,n_idx]
            if idx_prev != 0
                prev_ψ₁, prev_ψ₂ = @view ψ[:,idx_prev]
                tmp = g_tmp[mode] * sqrt_of_int[nₐ+1]
                dψ₁ += -tmp*prev_ψ₂
                dψ₂ += -tmp*prev_ψ₁
            end
            if idx_next != 0
                next_ψ₁, next_ψ₂ = @view ψ[:,idx_next]
                tmp = f_tmp[mode] * sqrt_of_int[nₐ+2]
                dψ₁ += -tmp*(next_ψ₂-L_av*next_ψ₁)
                dψ₂ += -tmp*(next_ψ₁-L_av*next_ψ₂)
            end
        end
        dψ[1,n_idx]=dψ₁-totalΓ*ψ₁
        dψ[2,n_idx]=dψ₂-totalΓ*ψ₂
    end
    return nothing
end