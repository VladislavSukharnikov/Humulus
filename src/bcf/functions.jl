# =============================================================================
# Elementary building blocks for nonstationary functions. 
# =============================================================================

"""
    phasecomb(t, ω, φ, u, v)

Return the phase combination

    u * exp(-im*x) - v * exp(im*x)

where

    x = ω*t - φ/2.

Equivalently,

    u*cis(-x) - v*cis(x)

The following parameter choices recover common functions:

| Function    | `u`       | `v`       |
|:------------|:----------|:----------|
| `cos(x)`    | `0.5`     | `-0.5`    |
| `sin(x)`    | `0.5im`   | `0.5im`   |
| `-sin(x)`   | `-0.5im`  | `-0.5im`  |
| `exp(im*x)` | `0.0`     | `-1.0`    |
"""
@inline function phasecomb(
                    t::Float64, 
                    ω::Float64, 
                    φ::Float64, 
                    u::ComplexF64, 
                    v::ComplexF64
                )::ComplexF64
                
    x = ω * t - φ/2
    z = cis(x)
    return u * conj(z) - v * z
end