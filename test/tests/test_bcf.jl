@testset "BCF" begin
    @info "BCF: predefined constructors..."

    BCF_RTOL = 1e-12
    BCF_ATOL = 1e-14

    test_cases = (
        ("1-mode random", random_bcf(100)),
        ("100-mode random", random_bcf(100)),
        ("One-mode squeezed", one_mode_squeezed_bcf(5.0, 1.0, 1.5, 0.0, 1.0)),
        ("Three-mode squeezed", three_mode_squeezed_bcf(5.0, 2.0, 1.0, 0.5, 0.0, 1.0)),
    )

    @testset "Predefined BCFs" begin
        @testset "$name" for (name, bcf) in test_cases
            N = n_modes(bcf)

            @testset "Type parameters" begin
                @test typeof(bcf).parameters[1] == N
                @test typeof(bcf).parameters[2] == typeof(bcf.f_vector)
                @test typeof(bcf).parameters[3] == typeof(bcf.g_vector)
            end

            @testset "Type stability" begin
                t, s = rand(2)

                @inferred bcf(t, s)

                for j in 1:N
                    @inferred bcf.f_vector[j](t)
                    @inferred bcf.f_vector[j](s)
                end
            end

            @testset "Fields" begin
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

            @testset "Return types" begin
                t, s = rand(2)

                @test bcf(t, s) isa ComplexF64

                for j in 1:N
                    @test bcf.f_vector[j](t) isa ComplexF64
                    @test bcf.g_vector[j](t) isa ComplexF64
                end
            end

            @testset "Physical properties" begin
                ts = 0:0.5:20

                @testset "Diagonal is real" begin
                    for t in ts
                        @test isapprox(imag(bcf(t, t)), 0.0; atol=BCF_ATOL)
                    end
                end

                @testset "Diagonal is positive" begin
                    for t in ts
                        @test real(bcf(t, t)) ≥ 0.0
                    end
                end

                @testset "Hermitian symmetry" begin
                    for i in eachindex(ts)
                        for j in 1:i
                            tᵢ, sⱼ = ts[i], ts[j]

                            @test isapprox(
                                bcf(tᵢ, sⱼ),
                                conj(bcf(sⱼ, tᵢ));
                                rtol=BCF_RTOL,
                                atol=BCF_ATOL,
                            )
                        end
                    end
                end
            end

            @testset "Zero allocations" begin
                t, s = rand(2)
                @test @ballocated($bcf($t, $s)) == 0
            end
        end
    end

    @info "BCF: construction..."

    @testset "Construction" begin
        for n in (2, 3, 4, 20, 30, 40)
            Γ = SVector{n,Float64}(rand() for _ in 1:n)
            ω = SVector{n,Float64}(rand() for _ in 1:n)
            G = SVector{n,Float64}(rand() for _ in 1:n)
            φ = SVector{n,Float64}(rand() for _ in 1:n)
            u = SVector{n,ComplexF64}(rand() + 1im * rand() for _ in 1:n)
            v = SVector{n,ComplexF64}(rand() + 1im * rand() for _ in 1:n)

            f_vector = SVector{n}(FuncWrapper(phasecomb, (ω[i], φ[i], u[i], v[i])) for i in 1:n)
            g_vector = SVector{n}(FuncWrapper(phasecomb, (ω[i], φ[i], u[i], v[i])) for i in 1:n)

            bcf = BCF{n,typeof(f_vector),typeof(g_vector)}(Γ, G, f_vector, g_vector)

            @test bcf.f_vector == f_vector
            @test bcf.g_vector == g_vector
            @test bcf.Γ == Γ
            @test bcf.G == G
            @test n_modes(bcf) == n
        end
    end

    @info "BCF: completed."
end;