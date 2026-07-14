# =============================================================================
# Function wrappers. 
# =============================================================================

"""
    FuncWrapper{F,P} <: Function

Callable object storing a function together with its parameters.

# Type parameters
- `F`: wrapped callable type.
- `P`: tuple type of the stored parameters.

# Fields
- `f`: wrapped callable.
- `params`: stored parameter tuple.
"""
struct FuncWrapper{F,P}<:Function
    f::F
    params::P
end

@inline (func::FuncWrapper)(t::Float64) = func.f(t, func.params...)


# =============================================================================
# Bath-correlation function. 
# =============================================================================

"""
    BCF{N,FVector,GVector} <: Function

Bath correlation function represented as a finite sum of exponentially
decaying modes,

```math
\\alpha(t,s) = \\sum_j G_j^2 f_j(t) g^*_j(s) e^{-\\Gamma_j |t-s|}.
```

# Type parameters
- `N`: number of modes.
- `FVector`: static vector of functions `fⱼ(t)`.
- `GVector`: static vector of functions `gⱼ(s)`.

# Fields
- `Γ`: decay rates.
- `G`: coupling amplitudes.
- `f_vector`: functions `fⱼ(t)`.
- `g_vector`: functions `gⱼ(s)`.
"""
struct BCF{N,
           FVector<:SVector{N,<:FuncWrapper}, 
           GVector<:SVector{N,<:FuncWrapper}
        } <: Function

    # Decay rates Γⱼ
    Γ :: SVector{N,Float64}

    # Coupling amplitudes Gⱼ
    G :: SVector{N,Float64}

    # Functions fⱼ(t)
    f_vector :: FVector

    # Functions gⱼ(s)
    g_vector :: GVector
end

# Number of effective modes.
n_modes(::BCF{N}) where {N} = N


# =============================================================================
# Bath-correlation function evaluation.
# =============================================================================

"""
    (bcf::BCF)(t, s) -> ComplexF64

Evaluate the bath-correlation function.

# Arguments
- `t`: first time.
- `s`: second time.

# Returns
- Value of the bath-correlation function at `(t, s)`.
"""
function (bcf::BCF{N})(
                    t::Float64, 
                    s::Float64,
                )::ComplexF64 where {N}
            
    (; Γ, G, f_vector, g_vector) = bcf

    # Time difference entering the exponential memory kernels.
    τ = abs(t-s)

    output::ComplexF64 = 0.0 + 0.0im

    for j in 1:N
        output += G[j]^2 * f_vector[j](t) * conj(g_vector[j](s)) * exp(-Γ[j] * τ)
    end

    return output
end