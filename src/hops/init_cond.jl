@inline function init_hops!(u0::ArrayPartition, c_g::T, c_e::T) where {T<:Number}
    ψ0, X0 = u0.x
    ψ0[1] = c_g # ψ_0[1,1]
    ψ0[2] = c_e # ψ_0[2,1]

    ψ0[3:end] .= 0.0im
    X0         .= 0.0im
    return nothing
end

function init_hops(fock_dim::Int, N::Int, c_g::T, c_e::T) where {T<:Number}

    ψ0 :: Array{ComplexF64,2} = zeros(ComplexF64, 2, fock_dim)
    X0 :: Array{ComplexF64,2} = zeros(ComplexF64, N, 1)

    u0 = ArrayPartition(ψ0, X0)
    init_hops!(u0, c_g, c_e)

    return u0
end