@testset "FockSpace" begin
   @info "FockSpace: construction..."

    key_types = (Int16, Int32, Int64)
    ind_types = (Int16, Int32, Int64)

    occupancy_cases = (
        ("scalar", 5),
        ("tuple", (3, 5)),
        ("vector", [3, 5]),
    )

    test_cases = (
        (
            "occupancy=$occupancy_title, key=$(nameof(KeyType)), index=$(nameof(IndType))",
            2,
            Val(2),
            max_occupancies,
            KeyType,
            IndType,
        )
        for (occupancy_title, max_occupancies) in occupancy_cases
        for KeyType in key_types
        for IndType in ind_types
    )

    @testset "$title" for (title, N, N_val, max_occupancies, KeyType, IndType) in test_cases
        fock_space = FockSpace(N_val, max_occupancies, KeyType, IndType)

        if max_occupancies isa Integer
            expected_dim = binomial(big(max_occupancies+N),N)
        else
            expected_dim = prod(max_occupancies.+1)
        end

        @testset "Type parameters" begin
            @test n_modes(fock_space) == N
            @test typeof(fock_space).parameters[1] == N
            @test typeof(fock_space).parameters[2] == KeyType
            @test typeof(fock_space).parameters[3] == IndType
        end

        @testset "Fields" begin
            if max_occupancies isa Vector
                @test fock_space.max_occupancies == Tuple(max_occupancies)
            else
                @test fock_space.fock_dim == expected_dim
            end
            @test fock_space.fock_dim == expected_dim

            @test eltype(fock_space.basis_states) == SVector{N,KeyType}
            @test fock_space.state_to_index isa Dict{SVector{N,KeyType},IndType}
            @test eltype(fock_space.raise_index) == IndType
            @test eltype(fock_space.lower_index) == IndType
        end

        @testset "Dimensions" begin
            @test length(fock_space.basis_states) == expected_dim
            @test length(fock_space.state_to_index) == expected_dim
            @test size(fock_space.raise_index) == (N, expected_dim)
            @test size(fock_space.lower_index) == (N, expected_dim)
        end

        @testset "State mapping" begin
            for (index, state) in pairs(fock_space.basis_states)
                @test fock_space.state_to_index[state] == index
            end
        end

        @testset "Occupation bounds" begin
            for state in fock_space.basis_states
                @test all(state .≥ 0)
                @test all(state .≤ max_occupancies)
            end
        end

        zero_index = zero(IndType)

        @testset "Creation indices" begin
            for (item, state) in enumerate(fock_space.basis_states)
                for mode in 1:N
                    raised_state = setindex(state, state[mode] + one(KeyType), mode)
                    raised_item  = fock_space.raise_index[mode, item]

                    if haskey(fock_space.state_to_index, raised_state)
                        expected_item = fock_space.state_to_index[raised_state]

                        @test raised_item == expected_item
                        @test raised_item != zero_index
                        @test fock_space.basis_states[raised_item] == raised_state
                    else
                        @test raised_item == zero_index
                    end
                end
            end
        end

        @testset "Annihilation indices" begin
            for (item, state) in enumerate(fock_space.basis_states)
                for mode in 1:N
                    lowered_item = fock_space.lower_index[mode, item]

                    if iszero(state[mode])
                        @test lowered_item == zero_index
                    else
                        lowered_state = setindex(state, state[mode] - one(KeyType), mode)

                        @test haskey(fock_space.state_to_index, lowered_state)

                        expected_item = fock_space.state_to_index[lowered_state]

                        @test lowered_item == expected_item
                        @test lowered_item != zero_index
                        @test fock_space.basis_states[lowered_item] == lowered_state
                    end
                end
            end
        end

        @testset "Inverse neighbors" begin
            for item in eachindex(fock_space.basis_states)
                for mode in 1:N
                    raised_item = fock_space.raise_index[mode, item]

                    if !iszero(raised_item)
                        @test fock_space.lower_index[mode, raised_item] == item
                    end

                    lowered_item = fock_space.lower_index[mode, item]

                    if !iszero(lowered_item)
                        @test fock_space.raise_index[mode, lowered_item] == item
                    end
                end
            end
        end
    end

    @testset "Invalid inputs" begin

        @test_throws MethodError FockSpace(Val(-2), (1, 5))
        @test_throws MethodError FockSpace(Val(0.1), (1, 5))
        @test_throws MethodError FockSpace(Val(0im), (1, 5))
        @test_throws MethodError FockSpace(Val(2), (5, 5, 5), Int16, Int32)
        @test_throws MethodError FockSpace(Val(2), [0.1, 0.2], Int16, Int32)

        @test_throws DomainError FockSpace(Val(0), 10)
        @test_throws DomainError FockSpace(Val(2), (-1, 5))
        @test_throws DomainError FockSpace(Val(2), [-1, 5])
        @test_throws DimensionMismatch FockSpace(Val(2), [5], Int16, Int32)
        @test_throws DimensionMismatch FockSpace(Val(2), [5,5,5], Int16, Int32)
        @test_throws DimensionMismatch FockSpace(Val(2), [5], Int16, Int32)
        @test_throws ArgumentError FockSpace(Val(1), 200, Int8, Int8)
        @test_throws ArgumentError FockSpace(Val(1), 200, Int, Int8)
        @test_throws ArgumentError FockSpace(Val(1), 200, Int8, Int)
    end

    @testset "Zero occupancy" begin
        fock_space = FockSpace(Val(2), (0, 0), Int16, Int16)

        @test fock_space.fock_dim == 1
        @test fock_space.basis_states == [SVector{2,Int16}(0, 0)]
        @test length(fock_space.state_to_index) == 1
        @test size(fock_space.raise_index) == (2, 1)
        @test size(fock_space.lower_index) == (2, 1)
    end

    @testset "Equivalent occupancy inputs" begin
        tuple_space  = FockSpace(Val(2), (5, 5), Int16, Int32)
        vector_space = FockSpace(Val(2), [5, 5], Int16, Int32)

        @test tuple_space.basis_states == vector_space.basis_states
        @test tuple_space.raise_index == vector_space.raise_index
        @test tuple_space.lower_index == vector_space.lower_index
    end

    @info "FockSpace: completed."
end;