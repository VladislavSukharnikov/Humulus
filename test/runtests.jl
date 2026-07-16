include("test_setup.jl")

@testset "Humulus" begin
    @info "Running Humulus tests..."

    include("tests/test_bcf.jl")
    include("tests/test_fock.jl")
    include("tests/test_hme.jl")
    include("tests/test_noise.jl")
    include("tests/test_hops.jl")

    @info "Finished Humulus tests."
end;