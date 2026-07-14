include("prof_utils.jl");

# =============================================================================
# Profile repeated FockSpace construction.
# =============================================================================

let 
    # Number of modes.
    N = 4;

    # Types used for the occupation-number lookup table.
    KeyType, IndType = Int16, Int32;

    # Maximum occupation number allowed for each mode.
    max_occupancies = (10, 10, 10, 10);

    # Trigger compilation.
    profile_fockspace(N, max_occupancies, KeyType, IndType; repetition_number=1);

    repetition_number = 1000;
    @profview profile_fockspace(N, max_occupancies, KeyType, IndType, repetition_number=repetition_number)
end