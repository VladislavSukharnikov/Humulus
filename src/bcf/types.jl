# =============================================================================
# Function wrappers. 
# =============================================================================

"""
    FuncWrapper{F,P} <: Function

Callable object that stores a function together with a tuple of parameters.

Calling `FuncWrapper(f, params)(t)` is equivalent to

    f(t, params...)

# Type parameters
- `F`: type of the wrapped callable.
- `P`: type of the stored parameter tuple.

# Fields
- `f`: wrapped callable.
- `params`: tuple of stored parameters.
"""
struct FuncWrapper{F,P}<:Function
    f      :: F
    params :: P
end

# Evaluate the wrapped function using the stored parameters.
@inline (func::FuncWrapper)(t::Float64) = func.f(t, func.params...)

# Value equality.
Base.:(==)(x::FuncWrapper, y::FuncWrapper) =
    x.f == y.f &&
    x.params == y.params

# Exact equality.
Base.isequal(x::FuncWrapper, y::FuncWrapper) =
    isequal(x.f, y.f) &&
    isequal(x.params, y.params)

    
# =============================================================================
# Bath-correlation function. 
# =============================================================================

"""
    BCF{N,FVector,GVector} <: Function

Callable bath-correlation function represented as a finite sum of exponentially
decaying modes,

```math
\\alpha(t,s) = \\sum_j G_j^2 f_j(t) g_j^*(s) e^{-\\Gamma_j |t-s|}.
```

# Type parameters
- `N`: number of modes.
- `FVector`: static vector of `FuncWrapper`s representing the functions `fⱼ(t)`.
- `GVector`: static vector of `FuncWrapper`s representing the functions `gⱼ(s)`.

# Fields
- `Γ`: mode decay rates.
- `G`: mode coupling amplitudes.
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

# Value equality.
Base.:(==)(x::BCF, y::BCF) =
    x.Γ == y.Γ &&
    x.G == y.G &&
    x.f_vector == y.f_vector &&
    x.g_vector == y.g_vector

# Exact equality.
Base.isequal(x::BCF, y::BCF) =
    isequal(x.Γ, y.Γ) &&
    isequal(x.G, y.G) &&
    isequal(x.f_vector, y.f_vector) &&
    isequal(x.g_vector, y.g_vector)


# =============================================================================
# Bath-correlation function evaluation.
# =============================================================================

"""
    (bcf::BCF)(t, s)

Evaluate the bath correlation function at the times `t` and `s`.

# Arguments
- `t::Float64`: first time.
- `s::Float64`: second time.
"""
function (bcf::BCF{N})(
                    t::Float64, 
                    s::Float64,
                )::ComplexF64 where {N}
            
    (; Γ, G, f_vector, g_vector) = bcf

    # Time difference entering the exponential memory kernels.
    τ = abs(t-s)

    output = 0.0im
    for j in 1:N
        output += G[j]^2 * f_vector[j](t) * conj(g_vector[j](s)) * exp(-Γ[j] * τ)
    end

    return output
end