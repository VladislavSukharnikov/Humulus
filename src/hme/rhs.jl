"""
    (hme::HME{N})(dρ, ρ, solver_params, t)

Evaluate the right-hand side of the hierarchy of master equations (HME).

The result is written in-place to `dρ`, overwriting its previous contents.

# Arguments
- `dρ`: output array containing the time derivative.
- `ρ`: hierarchy state.
- `solver_params`: solver parameters created by `create_solver_params`.
- `t::Float64`: time.
"""
function (hme::HME{N})(
                    dρ::Array{ComplexF64,4}, 
                    ρ::Array{ComplexF64,4}, 
                    solver_params::NamedTuple, 
                    t::Float64
                ) where {N}

    (; G, f_vector, g_vector) = solver_params
    (; f_tmp, g_tmp, _sqrt_of_int) = hme

    # Evaluate the time-dependent coefficients.
    for j in 1:N
        f_tmp[j] = 1im * G[j] * f_vector[j](t)
        g_tmp[j] = 1im * G[j] * conj(g_vector[j](t))
    end
 
    # Compute the HME right-hand side.
    _fill_dρ!(dρ, ρ, solver_params, f_tmp, g_tmp, _sqrt_of_int) 

    return nothing
end


"""
    _fill_dρ!(dρ, ρ, solver_params, f_tmp, g_tmp, sqrt_of_int)

Compute the right-hand side of the hierarchy of master equations (HME).

The result is written in-place to `dρ`, overwriting its previous contents.
"""
@inline function _fill_dρ!(
                        dρ::Array{ComplexF64,4},
                        ρ::Array{ComplexF64,4},
                        solver_params::NamedTuple,
                        f_tmp::MVector,
                        g_tmp::MVector,
                        sqrt_of_int::SVector
                    )
    (; N, ν₀, basis_states, fock_dim, Γ, raise_index, lower_index) = solver_params
    mode_range = 1:N

    # Iterate over all hierarchy elements.
    @inbounds for n_idx in 1:fock_dim
        for m_idx in n_idx:fock_dim

            # Current auxiliary density operator.
            rho_nm  = @view  ρ[:, :, n_idx, m_idx]
            drho_nm = @view dρ[:, :, n_idx, m_idx]

            n = basis_states[n_idx]
            m = basis_states[m_idx]

            ρ₁₁::ComplexF64 = rho_nm[1, 1]
            ρ₂₁::ComplexF64 = rho_nm[2, 1]
            ρ₁₂::ComplexF64 = rho_nm[1, 2]
            ρ₂₂::ComplexF64 = rho_nm[2, 2]

            dρ₁₁::ComplexF64 = 0.0im
            dρ₂₁::ComplexF64 =-1im * ν₀ * ρ₂₁
            dρ₁₂::ComplexF64 = 1im * ν₀ * ρ₁₂
            dρ₂₂::ComplexF64 = 0.0im

            totalΓ = 0.0

            # sqrt_of_int[k] = √(k-1).
            # Accumulate couplings to neighboring hierarchy elements
            for mode in mode_range
                nₐ = n[mode]
                mₐ = m[mode]
                totalΓ += Γ[mode] * (nₐ + mₐ)

                n_idx_prev, n_idx_next = lower_index[mode, n_idx], raise_index[mode, n_idx]
                m_idx_prev, m_idx_next = lower_index[mode, m_idx], raise_index[mode, m_idx]

                if n_idx_next !=0 
                    rho_n_next = @view ρ[:,:, n_idx_next, m_idx]
                    tmp = -f_tmp[mode]*sqrt_of_int[nₐ+2]

                    dρ₁₁ += tmp * (rho_n_next[2, 1] - rho_n_next[1, 2])
                    dρ₂₁ += tmp * (rho_n_next[1, 1] - rho_n_next[2, 2])
                    dρ₁₂ += tmp * (rho_n_next[2, 2] - rho_n_next[1, 1])
                    dρ₂₂ += tmp * (rho_n_next[1, 2] - rho_n_next[2, 1])
                end
                if n_idx_prev != 0
                    rho_n_prev = @view ρ[:, :, n_idx_prev, m_idx]
                    tmp = -g_tmp[mode] * sqrt_of_int[nₐ + 1]

                    dρ₁₁ += tmp * rho_n_prev[2, 1]
                    dρ₂₁ += tmp * rho_n_prev[1, 1]
                    dρ₁₂ += tmp * rho_n_prev[2, 2]
                    dρ₂₂ += tmp * rho_n_prev[1, 2]
                end
                if m_idx_next !=0 
                    rho_m_next = @view ρ[:, :, n_idx, m_idx_next]
                    tmp = -conj(f_tmp[mode]) * sqrt_of_int[mₐ + 2]

                    dρ₁₁ += tmp * (rho_m_next[1, 2] - rho_m_next[2, 1])
                    dρ₂₁ += tmp * (rho_m_next[2, 2] - rho_m_next[1, 1])
                    dρ₁₂ += tmp * (rho_m_next[1, 1] - rho_m_next[2, 2])
                    dρ₂₂ += tmp * (rho_m_next[2, 1] - rho_m_next[1, 2])
                end
                if m_idx_prev != 0
                    rho_m_prev = @view ρ[:, :, n_idx, m_idx_prev]
                    tmp = -conj(g_tmp[mode]) * sqrt_of_int[mₐ + 1]

                    dρ₁₁ += tmp * rho_m_prev[1, 2]
                    dρ₂₁ += tmp * rho_m_prev[2, 2]
                    dρ₁₂ += tmp * rho_m_prev[1, 1]
                    dρ₂₂ += tmp * rho_m_prev[2, 1]
                end
            end

            # Apply the damping contribution.
            drho_nm[1, 1] = dρ₁₁ - totalΓ * ρ₁₁
            drho_nm[2, 1] = dρ₂₁ - totalΓ * ρ₂₁
            drho_nm[1, 2] = dρ₁₂ - totalΓ * ρ₁₂
            drho_nm[2, 2] = dρ₂₂ - totalΓ * ρ₂₂

            if n_idx != m_idx
                dρ[1, 1, m_idx, n_idx] = conj(drho_nm[1, 1])
                dρ[2, 1, m_idx, n_idx] = conj(drho_nm[1, 2])
                dρ[1, 2, m_idx, n_idx] = conj(drho_nm[2, 1])
                dρ[2, 2, m_idx, n_idx] = conj(drho_nm[2, 2])
            end
            
        end
    end
    return nothing
end