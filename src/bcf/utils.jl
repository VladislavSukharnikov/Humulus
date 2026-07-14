# =============================================================================
# Testing utilities. 
# =============================================================================

"""
    create_random_bcf(N::Int)

Create a random nonstationary bath correlation function with `N` modes.

This helper is intended for tests only. The generated BCF uses
positive weights and is not meant to represent a specific physical system.
"""
function random_bcf(N::Int)
    @assert N > 0 "Number of modes should be positive and nonzero."

    if N > 100
        @warn "N=$N is large for StaticArrays; compilation may be slow."
    end

    return _random_bcf(Val(N))
end

# Internal implementation. `N` is passed through `Val` so that the
# StaticArrays size is known at compile time.
function _random_bcf(::Val{N}) where {N}
    Γ = SVector{N,Float64}(rand() for _ in 1:N)
    ω = SVector{N,Float64}(rand() for _ in 1:N)
    G = SVector{N,Float64}(rand() for _ in 1:N)

    φ = SVector{N,Float64}(rand() for _ in 1:N)
    u = SVector{N,ComplexF64}(rand(ComplexF64) for _ in 1:N)
    v = SVector{N,ComplexF64}(rand(ComplexF64) for _ in 1:N)

    f_vector = SVector{N}(FuncWrapper(phasecomb, (ω[i], φ[i], u[i], v[i])) for i in 1:N)
    g_vector = SVector{N}(FuncWrapper(phasecomb, (ω[i], φ[i], u[i], v[i])) for i in 1:N)

    return BCF{N, typeof(f_vector), typeof(g_vector)}(Γ, G, f_vector, g_vector)
end