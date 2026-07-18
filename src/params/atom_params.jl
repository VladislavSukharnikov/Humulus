# =============================================================================
# Atomic parameters. 
# =============================================================================

"""
    AtomParams

Parameters describing the two-level atom.

The initial atomic state is

```math
|ψ⟩ = c_g |g⟩ + c_e |e⟩.
```

# Fields
- `ν₀::Float64`: Atomic transition frequency.
- `c_g::ComplexF64`: Initial amplitude of the ground state.
- `c_e::ComplexF64`: Initial amplitude of the excited state.
"""
struct AtomParams
    ν₀::Float64
    c_g::ComplexF64
    c_e::ComplexF64
end