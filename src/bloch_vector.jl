"""
    bloch_vector(ρ_s)

Compute the expectation values of the Bloch-vector components from the
time-dependent reduced density matrix.

# Arguments
- `ρ_s::Array{ComplexF64,3}`: reduced density matrix sampled on the save-time grid.

# Returns
A tuple `(σˣ, σʸ, σᶻ)` containing the expectation values of the Pauli
operators at each saved time point.
"""
function bloch_vector(ρ_s::Array{ComplexF64,3})

    n_save = size(ρ_s)[end]

    σˣ = zeros(Float64,n_save)
    σʸ = zeros(Float64,n_save)
    σᶻ = zeros(Float64,n_save)

    σˣ .= 2 * real(ρ_s[1,2,:])
    σʸ .= 2 * imag(ρ_s[1,2,:])
    σᶻ .= real(ρ_s[2,2,:] .- ρ_s[1,1,:])

    return (σˣ, σʸ, σᶻ)
end