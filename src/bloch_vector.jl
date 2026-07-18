"""
    bloch_vector(ρ_s)

Return the expectation values of the Pauli operators for a time-dependent
reduced density matrix.

# Arguments
- `ρ_s::Array{ComplexF64,3}`: time-dependent reduced density matrix of size
  `(2, 2, n_save)`, where the third dimension corresponds to the save-time grid.
"""
function bloch_vector(ρ_s::Array{ComplexF64,3})

    n_save = size(ρ_s, 3)

    σˣ = zeros(Float64, n_save)
    σʸ = zeros(Float64, n_save)
    σᶻ = zeros(Float64, n_save)

    σˣ .= 2 * real(ρ_s[1,2,:])
    σʸ .= 2 * imag(ρ_s[1,2,:])
    σᶻ .= real(ρ_s[2,2,:] .- ρ_s[1,1,:])

    return (σˣ, σʸ, σᶻ)
end