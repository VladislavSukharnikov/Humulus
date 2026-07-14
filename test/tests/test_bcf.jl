testset_name = "BCF";
@testset "$testset_name" begin
    println("   $testset_name: setting configuration...")

    BCF_RTOL = 1e-12
    BCF_ATOL = 1e-14

    # Number of effective modes.
    N = 100
    bcf = random_bcf(N)


    println("   $testset_name: testing type parameters...")
    @testset "Type parameters" begin
        @test n_modes(bcf) == N
        @test typeof(bcf).parameters[1] == N
        @test typeof(bcf).parameters[2] == typeof(bcf.f_vector)
        @test typeof(bcf).parameters[3] == typeof(bcf.g_vector)
    end


    println("   $testset_name: testing type stability...")
    @testset "Type stability" begin
        t, s = rand(), rand()
        @inferred bcf(t, s)
    end


    println("   $testset_name: testing fields...")
    @testset "Expected fields" begin
        @test bcf.G isa SVector{N,Float64}
        @test bcf.Γ isa SVector{N,Float64}
        @test bcf.f_vector isa SVector{N,<:FuncWrapper}
        @test bcf.g_vector isa SVector{N,<:FuncWrapper}

        for j in 1:N
            @test bcf.f_vector[j] isa FuncWrapper
            @test bcf.g_vector[j] isa FuncWrapper
        end

        @test length(bcf.G) == N
        @test length(bcf.Γ) == N
        @test length(bcf.f_vector) == N
        @test length(bcf.g_vector) == N
    end
    

    println("   $testset_name: testing expected types...")
    @testset "Expected types" begin
        t, s = rand(), rand()
        @test bcf(t, s) isa ComplexF64

        for j in 1:N
            @test bcf.f_vector[j](t) isa ComplexF64
            @test bcf.g_vector[j](t) isa ComplexF64
        end
    end

    println("   $testset_name: physical validation...")
    @testset "Physical validation" begin
        ts = 0:0.5:50

        println("       $testset_name: real diagonal...")
        @testset "Real diagonal" begin
            for t in ts
                @test isapprox(imag(bcf(t, t)), 0.0; atol = BCF_ATOL)
            end
        end

        println("       $testset_name: positive diagonal...")
        @testset "Positive diagonal" begin
            for t in ts
                @test real(bcf(t, t)) ≥ 0.0
            end
        end

        println("       $testset_name: Hermitian symmetry...")
        @testset "Hermitian symmetry" begin
            for i in eachindex(ts)
                for j in 1:i
                    tᵢ = ts[i]
                    sⱼ = ts[j]

                    @test isapprox(
                        bcf(tᵢ, sⱼ),
                        conj(bcf(sⱼ, tᵢ));
                        rtol = BCF_RTOL,
                        atol = BCF_ATOL,
                    )
                end
            end
        end
    end


    println("   $testset_name: testing allocations...")
    t, s = rand(), rand()
    bench = @benchmark $bcf($t, $s)
    @test median(bench).memory==0


    println("   $testset_name: testing creation...")
    @testset "Testing creation" begin
        for n in [1,2,3,4,10,20,40]

            Γ = SVector{n,Float64}(rand() for _ in 1:n)
            ω = SVector{n,Float64}(rand() for _ in 1:n)
            G = SVector{n,Float64}(rand() for _ in 1:n)

            φ = SVector{n,Float64}(2*π*rand() for _ in 1:n)
            u = SVector{n,ComplexF64}(rand()+1im*rand() for _ in 1:n)
            v = SVector{n,ComplexF64}(rand()+1im*rand() for _ in 1:n)

            f_vector = SVector{n}(FuncWrapper(phasecomb, (ω[i], φ[i], u[i], v[i])) for i in 1:n)
            g_vector = SVector{n}(FuncWrapper(phasecomb, (ω[i], φ[i], u[i], v[i])) for i in 1:n)
                
            bcf_created = BCF{n, typeof(f_vector), typeof(g_vector)}(Γ, G, f_vector, g_vector)

            @test bcf_created.f_vector == f_vector
            @test bcf_created.g_vector == g_vector
            @test bcf_created.Γ == Γ
            @test bcf_created.G == G
            @test n_modes(bcf_created) == n
        end
    end


    println("   $testset_name: test complete.")
end