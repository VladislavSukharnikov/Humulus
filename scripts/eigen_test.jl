using Humulus
using LinearAlgebra
using Plots


grid_size = 10000
time_grid = Humulus.TimeGrid(0.0, 1.0, grid_size)
bcf_matrix = Matrix{ComplexF64}(undef, grid_size, grid_size)
    @time @inbounds for t in 1:grid_size
        for s in 1:t
            bcf_matrix[t, s] = bcf(time_grid[t], time_grid[s])
        end
    end

    t=rand(1:grid_size)
    s=rand(1:grid_size)
    bcf_matrix[t,s] 
    conj(bcf_matrix[s,t]) 


    bcf_matrix1 = Hermitian(bcf_matrix, :L)


    @time out = eigen!(bcf_matrix1);

    @time out = cholesky!(bcf_matrix1);


# # # ## # 


bcf = let
    ω₀ = 5.0        # central (carrier) frequency
    Γ₀ = 2.0        # input spectral width
    κ  = 1.0        # cavity loss rate
    ϵ  = 0.5        # effective DPA pump amplitude
    φ  = Float64(π) # squeezing phase
    γ  = 1.0        # atom–field coupling strength

    bcf = three_mode_squeezed_bcf(ω₀, Γ₀, κ, ϵ, φ, γ)
end

t_end = 1.0;
grid_size = 4000;

# Construct the discretization grid.
time_grid = Humulus.TimeGrid(0.0, t_end, grid_size)

# Assemble the discretized BCF kernel matrix.
bcf_matrix = Matrix{ComplexF64}(undef, grid_size, grid_size)
@inbounds for t in 1:grid_size, s in 1:grid_size
    bcf_matrix[t,s] = bcf(time_grid[t], time_grid[s])
end

# Compute the eigendecomposition of the BCF kernel.
@time vals, vecs = eigen(Hermitian(bcf_matrix));
@time vals, vecs = eigen!(Hermitian(bcf_matrix));


# # # # # 
 # MKL

# Assemble the discretized BCF kernel matrix.
bcf_matrix = Matrix{ComplexF64}(undef, grid_size, grid_size)
@inbounds for t in 1:grid_size, s in 1:grid_size
    bcf_matrix[t,s] = bcf(time_grid[t], time_grid[s])
end

# Compute the eigendecomposition of the BCF kernel.
@time vals, vecs = eigen(Hermitian(bcf_matrix));
@time vals, vecs = eigen!(Hermitian(bcf_matrix));
