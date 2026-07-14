include("bench_utils.jl");

# =============================================================================
# Benchmark FockSpace construction. 
# =============================================================================

let
    title = "Benchmark FockSpace construction." 
    print_header(title)

    N = 3

    KeyType = Int32
    IndType = Int16

    max_occupancies = (10, 10, 10)

    @info "Benchmark parameters:" N KeyType IndType max_occupancies

    fock_space = Humulus.FockSpace(Val(N), max_occupancies, KeyType, IndType)

    bench = @benchmark Humulus.FockSpace($(Val(N)), $max_occupancies, $KeyType, $IndType)
    benchmark_construction(fock_space, bench)
end