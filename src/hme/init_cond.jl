"""
    init_hme!(ρ, c_g, c_e)

Initialize the HME hierarchy in place.

The physical density matrix is initialized to the pure atomic state

```math
|ψ⟩ = c_g |g⟩ + c_e |e⟩,
```

while all auxiliary density operators are initialized to zero.

# Arguments
- `ρ`: hierarchy state.
- `c_g`: ground-state amplitude.
- `c_e`: excited-state amplitude.
"""
function init_hme!(ρ::AbstractArray{ComplexF64,4}, c_g::T, c_e::T) where {T<:Number}
    # Clear all entries before writing the initial state.
    fill!(ρ, 0.0im)

    ρ[1,1,1,1] = abs2(c_g)
    ρ[2,1,1,1] = c_e * conj(c_g)
    ρ[1,2,1,1] = c_g * conj(c_e)
    ρ[2,2,1,1] = abs2(c_e)

    return nothing
end


"""
    init_hme(fock_dim, c_g, c_e)

Allocate and initialize the HME hierarchy.

The physical density matrix is initialized to the pure atomic state
specified by `c_g` and `c_e`, while all auxiliary density operators are
initialized to zero.

# Arguments
- `fock_dim`: dimension of the pseudo-Fock space.
- `c_g`: ground-state amplitude.
- `c_e`: excited-state amplitude.
"""
function init_hme(fock_dim::Int, c_g::T, c_e::T) where {T<:Number}
    # Allocate the density matrix and initialize it in place.
    u0 = zeros(ComplexF64, 2, 2, fock_dim, fock_dim)
    init_hme!(u0, c_g, c_e)

    return u0
end