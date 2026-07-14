# =============================================================================
# Atomic parameters. 
# =============================================================================

"""
    AtomParams

Atomic parameters.

The initial atomic state is

```math
|ψ⟩ = c_g |g⟩ + c_e |e⟩.
```

# Fields
- `ν₀::Float64`: atomic transition frequency.
- `c_g::ComplexF64`: ground-state amplitude.
- `c_e::ComplexF64`: excited-state amplitude.
"""
struct AtomParams
    ν₀::Float64
    c_g::ComplexF64
    c_e::ComplexF64
end