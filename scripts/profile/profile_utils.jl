using Humulus

# =============================================================================
# Utilities.
# =============================================================================


"""
    profile_bcf_eval(bcf::BCF, t::Float64, s::Float64, repetition_number::Int=1)

Evaluates `bcf(t, s)` repeatedly and accumulates the result. This is
intended for profiling, where a single BCF evaluation is too fast for
Julia's sampling profiler to collect meaningful data.
"""
function profile_bcf_eval(bcf::Humulus.BCF, t::Float64, s::Float64; repetition_number::Int=1)
    # Accumulator for the repeated BCF evaluations.
    out = 0.0im

    # Repeat many times so that profiling has enough runtime to sample.
    for _ in 1:repetition_number
        out+=bcf(t,s)
    end

    # Return the accumulated value so the loop cannot be trivially ignored.
    return out
end


"""
    profile_fockspace(repetition_number::Int, N::Int, max_occupancies::Union{Int, Tuple},
                           KeyType::Type{<:Integer}, IndType::Type{<:Integer})

Constructs a `FockSpace` repeatedly. This is intended for profiling,
where constructing a single `FockSpace` is too fast for Julia's sampling
profiler to collect meaningful data.
"""
function profile_fockspace(
                        N::Int, 
                        max_occupancies::Union{Int,Tuple}, 
                        KeyType::Type{<:Integer}, 
                        IndType::Type{<:Integer}; 
                        repetition_number::Int=1
                    )
    # Accumulator for the repeated construction.
    out = 0

    # Repeat many times so that profiling has enough runtime to sample.
    for _ in 1:repetition_number
        fock_space = Humulus.FockSpace(Val(N), max_occupancies, KeyType, IndType)
        out+=fock_space.fock_dim
    end

    # Return the accumulated value so the loop cannot be trivially ignored.
    return out
end


"""
    profile_hme_eval!(hme!::HME, dρ::Array{ComplexF64,4}, ρ::Array{ComplexF64,4},
                      solver_params::NamedTuple, t::Float64,
                      repetition_number::Int=1)

Evaluates `hme!` repeatedly. This is intended for profiling, where a
single call to `hme!` is too fast for Julia's sampling profiler to collect
meaningful data.

Returns the number of completed evaluations.
"""
function profile_hme_eval!(
                        hme!::Humulus.HME, 
                        dρ::Array{ComplexF64,4}, 
                        ρ::Array{ComplexF64,4}, 
                        solver_params::NamedTuple, 
                        t::Float64; 
                        repetition_number::Int=1
                    )
                    
    # Accumulator for the repeated hme!(...) evaluations.
    out = 0

    # Repeat many times so that profiling has enough runtime to sample.
    for _ in 1:repetition_number
        out+=1
        hme!(dρ, ρ, solver_params, t)
    end

    # Return the accumulated value so the loop cannot be trivially ignored.
    return out
end