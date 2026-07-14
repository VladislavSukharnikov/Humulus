include("debug_utils.jl");

# =============================================================================
# Debugging: @code_warntype inspection of one_mode_squeezed_bcf(...).
# =============================================================================

let 
    # Construct BCF.
    bcf = one_mode_squeezed_bcf(5.0, 1.0, 1.5, 0.0, 1.0)
    t,s = rand(2)

    # Trigger compilation.
    bcf(t,s)

    @code_warntype bcf(t,s)
end


# =============================================================================
# Debugging: @code_warntype inspection of three_mode_squeezed_bcf(...).
# =============================================================================
let 
    # Construct BCF.
    bcf = three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)
    t,s = rand(2)

    # Trigger compilation.
    bcf(t,s)

    @code_warntype bcf(t,s)
end