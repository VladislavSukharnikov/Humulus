testset_name = "FockSpace";
@testset "$testset_name" begin
    println("   $testset_name: setting configuration...")

    # Number of effective modes.
    N = 4

    KeyType = Int16
    IndType = Int32

    # Maximum occupation number allowed for each mode.
    max_occupancies = (9, 9, 9, 9)

    fock_space = FockSpace(Val(N), max_occupancies, KeyType, IndType)
    @test fock_space.max_occupancies == max_occupancies

    println("   $testset_name: testing type parameters...")
    @testset "Type parameters" begin
        @test n_modes(fock_space) == N
        @test typeof(fock_space).parameters[1] == N
        @test typeof(fock_space).parameters[2] == KeyType
        @test typeof(fock_space).parameters[3] == IndType
    end

    println("   $testset_name: testing type stability...")
    @testset "Type stability" begin
        @test (@inferred FockSpace(Val(N), max_occupancies, KeyType, IndType)) isa FockSpace
    end

    println("   $testset_name: testing dimensionality...")
    @testset "Dimensionality" begin
        @test length(fock_space.basis_states) == 10^4
        @test fock_space.fock_dim == 10^4
        @test size(fock_space.raise_index) == (4, 10^4)
        @test size(fock_space.lower_index) == (4, 10^4)
    end

    println("   $testset_name: testing fields...")
    @testset "Expected fields" begin
        for n in 1:4
            for max_occupancy in 1:10
                fock_space_generated = FockSpace(Val(n), max_occupancy, KeyType, IndType)
                expected_dim = binomial(big(max_occupancy+n),n)

                @test n_modes(fock_space_generated) == n
                @test fock_space_generated.fock_dim == expected_dim
                @test size(fock_space_generated.raise_index) == (n, expected_dim)
                @test size(fock_space_generated.lower_index) == (n, expected_dim)
                @test length(fock_space_generated.basis_states) == expected_dim
                @test length(fock_space_generated.state_to_index) == expected_dim
                @test length(fock_space_generated.state_to_index) == expected_dim
                @test eltype(fock_space_generated.basis_states) == SVector{n, KeyType}
                @test fock_space_generated.state_to_index isa Dict{SVector{n, KeyType}, IndType} 
            end
            max_occupancies_rand = NTuple{n,Int}(rand(1:4,n))
            fock_space_generated = FockSpace(Val(n), max_occupancies_rand, KeyType, IndType)
            expected_dim = prod(max_occupancies_rand.+1)

            @test n_modes(fock_space_generated) == n
            @test fock_space_generated.fock_dim == expected_dim
            @test size(fock_space_generated.raise_index) == (n, expected_dim)
            @test size(fock_space_generated.lower_index) == (n, expected_dim)
            @test length(fock_space_generated.basis_states) == expected_dim
            @test length(fock_space_generated.state_to_index) == expected_dim
            @test length(fock_space_generated.state_to_index) == expected_dim
            @test eltype(fock_space_generated.basis_states) == SVector{n, KeyType}
            @test fock_space_generated.state_to_index isa Dict{SVector{n, KeyType}, IndType} 
        end
    end
    println("   $testset_name: test complete.")
end