"""
    sampler_from_cache(path::String; checks=true) -> NoiseSampler

Load precomputed noise data from `path`.

The file must contain a serialized `BCFEigen` object. When `checks=true`
(the default), the loaded data are validated before constructing the
returned `NoiseSampler`.

# Arguments
- `path`: path to a serialized `BCFEigen` object.

# Keyword arguments
- `checks`: whether to validate the loaded data before constructing the
  noise generator.

# Returns
- `NoiseSampler`: noise generator initialized from the stored time grid,
  eigenvectors, and square roots of the BCF kernel eigenvalues.

# Throws
- An exception if the file cannot be read.
- `AssertionError` if `checks=true` and the loaded data are invalid.
"""
function sampler_from_cache(
                path::String; 
                checks::Bool=true, 
                clear_cache::Bool=false
            )

    @assert isfile(path) "Cache file does not exist: $path"

    # Load precomputed noise parameters from disk.
    bcf_eigen::BCFEigen = load_object(path)::BCFEigen
    if clear_cache
        rm(path)
    end

    # Extract the stored eigendecomposition of the BCF.
    (; time_grid, vals, vecs) = bcf_eigen
    n_points = time_grid.n_points
    if checks
        @assert length(vals) == n_points "Expected $n_points eigenvalues."
        @assert size(vecs) == (n_points, n_points) "Expected a $n_points × $n_points eigenvector matrix."
    end
        
    # Warn if any eigenvalues are negative. This may indicate numerical
    # inaccuracies or that the BCF is not positive semidefinite.
    if any(vals .< 0)
        @warn "The discretized BCF matrix has negative eigenvalues."
    end

    # Workspace used during noise generation to avoid repeated allocations.
    container = zeros(ComplexF64, n_points)

    return NoiseSampler(time_grid, vals, vecs, container)
end


"""
    create_noise_cache(bcf, t_end, grid_size) -> String

Create or reuse a cached `BCFEigen` object.

The cache is identified by the bath correlation function, the final time,
and the discretization grid size. If a matching cache file already exists,
its path is returned. Otherwise, a new `BCFEigen` object is constructed,
saved to disk, and its path is returned.

Cache files are stored in the `noise_cache` directory.

# Arguments
- `bcf`: bath correlation function.
- `t_end`: final time of the discretization interval.
- `grid_size`: number of discretization points.

# Returns
- `String`: path to the cache file containing the corresponding
  `BCFEigen` object.

# Notes
An existing cache file is reused only if its contents are consistent with
the requested parameters. Otherwise, it is replaced with a newly generated
cache.
"""
function create_noise_cache(
                        bcf::BCF,
                        t_end::Float64,
                        grid_size::Int,
                    )

    dir = "noise_cache"
    mkpath(dir)

    key      = (t_end, bcf, grid_size)
    filename = string(hash(key), ".jld2")
    path     = joinpath(dir, filename)

    if isfile(path)
        try
            bcf_eigen::BCFEigen = load_object(path)::BCFEigen
            loaded_bcf::BCF     = bcf_eigen.bcf::BCF

            @info "Found existing cache."
            if  loaded_bcf == bcf &&
                bcf_eigen.time_grid[1] == 0.0 &&
                bcf_eigen.time_grid[end] >= t_end &&
                bcf_eigen.time_grid.n_points == grid_size 

                @info "Reusing noise cache."
                return path
            end
            @info "Creating noise cache."
        catch err
            @info "Failed to read noise cache. Recreating." exception=(err, catch_backtrace())
        end
    else
        @info "Rebuilding noise cache."
    end

    bcf_eigen = BCFEigen(bcf, t_end, grid_size)
    save_object(path, bcf_eigen)

    return path
end


function clear_noise_cache()
    dir = "noise_cache"
    isdir(dir) || return

    for entry in readdir(dir; join=true)
        rm(entry; recursive=true, force=true)
    end
end