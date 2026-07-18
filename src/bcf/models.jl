# =============================================================================
# Constructors of BCFs for physical models. 
# =============================================================================

"""
    three_mode_squeezed_bcf(ω₀, Γ₀, κ, ϵ, φ, γ)

Construct and return a three-mode effective bath-correlation function (BCF)
describing the output field of a degenerate parametric amplifier (DPA).

# Arguments
- `ω₀::Float64`: central output frequency.
- `Γ₀::Float64`: input linewidth.
- `κ::Float64`: cavity loss rate.
- `ϵ::Float64`: effective pump amplitude.
- `φ::Float64`: squeezing phase.
- `γ::Float64`: system–field coupling strength.

# Physical constraints

The parameters must satisfy

- `Γ₀ ≥ 0`, `κ ≥ 0`, `ϵ ≥ 0`, and `γ ≥ 0`;
- `κ > ϵ` (stable DPA);
- `Γ₀ > κ + ϵ` (the input bandwidth exceeds all effective linewidths).

If any of these conditions are violated, an exception is thrown.
"""
function three_mode_squeezed_bcf(
                            ω₀::Float64, 
                            Γ₀::Float64, 
                            κ::Float64, 
                            ϵ::Float64, 
                            φ::Float64, 
                            γ::Float64
                        )

    # Sanity checks for physical validity.

    Γ₀ >= 0 ||
        throw(DomainError(Γ₀, "The input linewidth Γ₀ must be non-negative."))

    κ >= 0 ||
        throw(DomainError(κ, "The cavity loss rate κ must be non-negative."))

    γ >= 0 ||
        throw(DomainError(γ, "The light–matter coupling γ must be non-negative."))

    ϵ >= 0 ||
        throw(DomainError(ϵ, "The pump amplitude ϵ must be non-negative. Use φ to specify the pump phase."))

    κ > ϵ ||
        throw(ArgumentError("The stability condition κ > ϵ must hold (got κ = $κ, ϵ = $ϵ)."))

    Γ₀ > κ + ϵ ||
        throw(ArgumentError("The input linewidth must satisfy Γ₀ > κ + ϵ (got Γ₀ = $Γ₀, κ + ϵ = $(κ + ϵ))."))

    # Effective mode decay rates.
    Γ₋ = κ - ϵ
    Γ₊ = κ + ϵ
    Γ  = SVector{3,Float64}(Γ₀, Γ₋, Γ₊)

    # Effective system–mode coupling strengths.
    G₀ = sqrt(γ * Γ₀/2)
    G₋ = sqrt(2 * γ * κ * ϵ/Γ₋ * Γ₀^2/(Γ₀^2 - Γ₋^2))
    G₊ = sqrt(2 * γ * κ * ϵ/Γ₊ * Γ₀^2/(Γ₀^2 - Γ₊^2))
    G  = SVector{3,Float64}(G₀, G₋, G₊)

    # Bogoliubov coefficients.
    u::ComplexF64 = (Γ₀^2 - κ^2 - ϵ^2)/sqrt((Γ₀^2 - Γ₊^2) * (Γ₀^2 - Γ₋^2))
    v::ComplexF64 = 2 * κ * ϵ/sqrt((Γ₀^2 - Γ₊^2)*(Γ₀^2 - Γ₋^2))

    # Time-dependent nonstationary functions.
    f_vector = SVector{3}(
        FuncWrapper(phasecomb, (ω₀, φ, u, v)),
        FuncWrapper(phasecomb, (ω₀, φ, 0.5 + 0im, -0.5 + 0im)), # cos(ω₀t+φ)
        FuncWrapper(phasecomb, (ω₀, φ, 0.5im, 0.5im)),          # sin(ω₀t+φ)
    )

    g_vector = SVector{3}(
        FuncWrapper(phasecomb, (ω₀, φ, u, v)),
        FuncWrapper(phasecomb, (ω₀, φ, 0.5 + 0im, -0.5 + 0im)), # cos(ω₀t+φ)
        FuncWrapper(phasecomb, (ω₀, φ, -0.5im, -0.5im)),        # -sin(ω₀t+φ)
    )

    return BCF{3, typeof(f_vector), typeof(g_vector)}(Γ, G, f_vector, g_vector)
end


"""
    one_mode_squeezed_bcf(ω₀, Γ, r, φ, γ)

Construct and return the bath-correlation function (BCF) for a single-mode
squeezed reservoir.

The resulting BCF consists of a single pseudo-mode with central frequency `ω₀`,
spectral half-width `Γ`, and effective coupling strength determined by `γ`. The
squeezing is characterized by the parameter `r` and phase `φ`.

# Arguments
- `ω₀::Float64`: central reservoir frequency.
- `Γ::Float64`: spectral half-width (memory decay rate).
- `r::Float64`: squeezing parameter.
- `φ::Float64`: squeezing phase.
- `γ::Float64`: system–field coupling strength.

# Physical constraints

The parameters must satisfy

- `Γ ≥ 0`, `r ≥ 0`, and `γ ≥ 0`.

If this condition is violated, an exception is thrown.
"""
function one_mode_squeezed_bcf(
                        ω₀::Float64, 
                        Γ::Float64, 
                        r::Float64, 
                        φ::Float64, 
                        γ::Float64,
                    )

    Γ >= 0 ||
        throw(DomainError(Γ, "The memory decay rate Γ must be non-negative."))

    γ >= 0 ||
        throw(DomainError(γ, "The light-matter coupling γ must be non-negative."))

    r >= 0 ||
        throw(DomainError(r, "The squeezing parameter r must be non-negative. Use φ to specify the squeezing phase."))

    # Effective coupling strength.
    Γ′ = SVector{1}(Γ)
    G  = SVector{1,Float64}(sqrt(γ*Γ/2))

    # Bogoliubov coefficients describing the squeezing transformation.
    u::ComplexF64 = cosh(r)
    v::ComplexF64 = sinh(r)

    # Time-dependent coefficients appearing in the BCF.
    f_vector = SVector{1}(FuncWrapper(phasecomb, (ω₀, φ, u, v)))
    g_vector = SVector{1}(FuncWrapper(phasecomb, (ω₀, φ, u, v)))

    return BCF{1,typeof(f_vector),typeof(g_vector)}(Γ′, G, f_vector, g_vector)
end