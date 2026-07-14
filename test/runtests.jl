using Test

@testset "Humulus" begin
    println("Running Humulus tests...")

    include("setup.jl")

    include("tests/test_bcf.jl")
    include("tests/test_fock_space.jl")
    include("tests/test_hme.jl")
    include("tests/test_noise.jl")
    include("tests/test_hops.jl")

    println("Finished Humulus tests.")
end;