"""
    init_hops!(u0, c_g, c_e)

Initialize the HOPS state in place.

The physical state is initialized to the pure atomic state

```math
|ψ⟩ = c_g |g⟩ + c_e |e⟩,
```

while all auxiliary hierarchy states and memory variables are initialized
to zero.

# Arguments
- `u0`: HOPS state.
- `c_g`: ground-state amplitude.
- `c_e`: excited-state amplitude.
"""
@inline function init_hops!(u0::ArrayPartition, c_g::T, c_e::T) where {T<:Number}
    ψ0, X0 = u0.x

    # Initialize the physical state.
    ψ0[1] = c_g
    ψ0[2] = c_e

    # Initialize the auxiliary hierarchy states and memory variables.
    ψ0[3:end] .= 0.0im
    X0 .= 0.0im

    return nothing
end

"""
    init_hops(fock_dim, N, c_g, c_e)

Allocate and initialize the HOPS state.

The physical state is initialized to the pure atomic state specified by
`c_g` and `c_e`, while all auxiliary hierarchy states and memory variables
are initialized to zero.

# Arguments
- `fock_dim`: dimension of the pseudo-Fock space.
- `N`: number of effective modes.
- `c_g`: ground-state amplitude.
- `c_e`: excited-state amplitude.
"""
function init_hops(fock_dim::Int, N::Int, c_g::T, c_e::T) where {T<:Number}

    ψ0::Array{ComplexF64,2} = zeros(ComplexF64, 2, fock_dim)
    X0::Array{ComplexF64,2} = zeros(ComplexF64, N, 1)

    # Allocate the HOPS state and initialize it in place.
    u0 = ArrayPartition(ψ0, X0)
    init_hops!(u0, c_g, c_e)

    return u0
end