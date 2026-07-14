include("prof_utils.jl")

# =============================================================================
# Profile repeated evaluation of a random BCF.
# =============================================================================

let 
    # Number of modes.
    N = 100

    # Construct a random BCF.
    bcf = Humulus.random_bcf(N)

    # Initialize two random time points.
    t, s = rand(2)

    # Trigger compilation.
    profile_bcf_eval(bcf, t, s; repetition_number=1)

    # Get profile.
    repetition_number = 1_000_000;
    @profview profile_bcf_eval(bcf, t, s; repetition_number)
end


# =============================================================================
# Profile repeated BCF evaluation for the single-mode squeezed reservoir.
# =============================================================================

let 
    # Construct BCF.
    bcf = one_mode_squeezed_bcf(5.0, 1.0, 1.5, 0.0, 1.0)

    # Initialize two random time points.
    t, s = rand(2)

    # Trigger compilation.
    profile_bcf_eval(bcf, t, s)

    # Get profile.
    repetition_number = 1_000_000
    @profview profile_bcf_eval(bcf, t, s; repetition_number=repetition_number)
end


# =============================================================================
# Profile repeated BCF evaluation for the three-mode squeezed reservoir.
# =============================================================================

let 
    # Construct BCF.
    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, Float64(π), 1.0)

    # Initialize two random time points.
    t, s = rand(), rand()

    # Trigger compilation.
    profile_bcf_eval(bcf, t, s; repetition_number=1);

    # Get profile.
    repetition_number = 1_000_000;
    @profview profile_bcf_eval(bcf, t, s, repetition_number=repetition_number)
end