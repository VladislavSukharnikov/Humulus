# =============================================================================
# Constructors of BCFs for physical models. 
# =============================================================================

"""
    three_mode_squeezed_bcf(ω₀, Γ₀, κ, ϵ, φ, γ) -> BCF

Construct the three-mode effective bath correlation function (BCF) describing the
output field of a degenerate parametric amplifier (DPA).

# Arguments
- `ω₀::Float64`: output central frequency.
- `Γ₀::Float64`: input Lorentzian width.
- `κ::Float64`: cavity loss rate.
- `ϵ::Float64`: effective DPA pump amplitude.
- `φ::Float64`: effective squeezing phase.
- `γ::Float64`: system–field coupling strength.

# Physical constraints
The parameters must satisfy
- `κ > ϵ` (stable DPA),
- `Γ₀ > κ + ϵ` (the input bandwidth is the largest linewidth),
- all rates and couplings are non-negative.
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

    @assert Γ₀ ≥ 0 "Input linewidth Γ₀ must be non-negative."
    @assert κ ≥ 0 "Cavity loss rate κ must be non-negative."
    @assert γ ≥ 0 "Light–matter coupling γ must be non-negative."
    @assert ϵ ≥ 0 "Pump amplitude ϵ must be non-negative (its phase is specified by φ)."

    @assert κ > ϵ "The stability condition κ > ϵ must hold."
    @assert Γ₀ > κ + ϵ "The input linewidth must satisfy Γ₀ > κ + ϵ."

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
    one_mode_squeezed_bcf(ω₀, Γ, r, φ, γ) -> BCF

Construct the bath correlation function (BCF) for a single-mode squeezed reservoir.

The resulting BCF consists of a single pseudo-mode with central frequency
`ω₀`, spectral half-width `Γ`, and effective coupling strength determined by
`γ`. The squeezing is characterized by the parameter `r` and phase `φ`.

# Arguments
- `ω₀::Float64`: central frequency of the reservoir.
- `Γ::Float64`: spectral half-width (memory decay rate).
- `r::Float64`: squeezing parameter.
- `φ::Float64`: squeezing phase.
- `γ::Float64`: atom–field coupling strength.

# Returns
- `BCF`: bath correlation function representing a one-mode squeezed reservoir.
"""
function one_mode_squeezed_bcf(
                        ω₀::Float64, 
                        Γ::Float64, 
                        r::Float64, 
                        φ::Float64, 
                        γ::Float64,
                    )

    @assert Γ>=0  "Memory decay rate `Γ` cannot be negative."
    @assert γ>=0  "Light-matter coupling `γ` cannot be negative."
    @assert r>=0  "Squeezing parameter `r` cannot be negative (its phase is specified by φ)."

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