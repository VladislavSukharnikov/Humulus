"""
    sampler_from_cache(path; checks=true, clear_cache=false)

Load a precomputed bath-correlation-function (BCF) eigendecomposition from
`path` and construct a `NoiseSampler`.

The file must contain a serialized `BCFEigen` object. When `checks=true`
(the default), the loaded data are validated before constructing the
returned sampler.

# Arguments
- `path::String`: path to a serialized `BCFEigen` object.

# Keyword arguments
- `checks::Bool`: whether to validate the loaded data.
- `clear_cache::Bool`: whether to remove the cache file after loading it.

# Exceptions

An exception is thrown if

- the cache file cannot be read;
- `checks=true` and the loaded data are invalid.
"""
function sampler_from_cache(
                        path::String;
                        checks::Bool=true,
                        clear_cache::Bool=false,
                    )

    @assert isfile(path) "Cache file does not exist: $path"

    cache = load_object(path)

    if checks && !(cache isa AbstractBCFDecomposition)
        throw(ArgumentError("Expected a subtype of AbstractBCFDecomposition, got $(typeof(cache))."))
    end

    clear_cache && rm(path)

    factor = _noise_factor(cache, checks)

    container = zeros(ComplexF64, cache.time_grid.n_points)

    return NoiseSampler(cache.time_grid, factor, container)
end

function _noise_factor(cache::BCFEigen, checks::Bool)
    n_points = cache.time_grid.n_points
    checks && @assert size(cache.vecs) == (n_points, n_points) "Expected a $n_points × $n_points matrix."
    return cache.vecs
end

function _noise_factor(cache::BCFCholesky, checks::Bool)
    n_points = cache.time_grid.n_points
    checks && @assert size(cache.chol) == (n_points, n_points) "Expected a $n_points × $n_points matrix."
    return cache.chol
end


"""
    get_covariance_cache(bcf, t_end, grid_size)

Return the path to a cached `BCFEigen` object.

The cache is identified by the bath correlation function, the final time,
and the discretization grid size. If a matching cache file already exists,
its path is returned. Otherwise, a new `BCFEigen` object is constructed,
saved to disk, and its path is returned.

Cache files are stored in the `bcf_eigen_cache` directory.

# Arguments
- `bcf::BCF`: bath correlation function.
- `t_end::Float64`: final time of the discretization interval.
- `grid_size::Int`: number of discretization points.

# Keyword arguments

- `logging::Bool=true`: if `true` (default), emit informational messages
  describing cache reuse, cache creation, and eigendecomposition progress.

# Notes

An existing cache file is reused only if its contents are consistent with
the requested parameters. Otherwise, it is replaced with a newly generated
cache.
"""
function get_covariance_cache(
                        ::Type{T},
                        bcf::BCF,
                        t_end::Float64,
                        grid_size::Int;
                        logging::Bool=true,
                    ) where {T<:AbstractBCFDecomposition}
    # Input validation.
    grid_size >= 2 ||
        throw(ArgumentError("`grid_size` must be at least 2, got $grid_size."))

    t_end > 0.0 ||
        throw(ArgumentError("`t_end` must be positive, got $t_end."))

    # Constructing cache filename.
    dir = "covariance_cache"
    mkpath(dir)
    filename = string(nameof(T), "_", bcf_hash(bcf), "_", hash((grid_size, t_end)), ".jld2")
    path     = joinpath(dir, filename)

    # Searching if the file already exists.
    if isfile(path)
        try
            cached::T = load_object(path)::T
            loaded_bcf::BCF = cached.bcf::BCF

            logging && @info "Found an existing cache."
            if loaded_bcf == bcf &&
                cached.time_grid[1] == 0.0 &&
                cached.time_grid[end] >= t_end &&
                cached.time_grid.t_end/cached.time_grid.n_points <= t_end/grid_size

                logging && @info "The existing cache is compatible. Reusing it."
                return path
            end
            logging && @info "The existing cache is incompatible. Rebuilding."
        catch err
            logging && @info "Failed to read the existing cache. Rebuilding." exception=(err, catch_backtrace())
        end
    else
        logging && @info "No existing cache found. Building a new one."
    end

    logging && @info "Computing the decomposition of the covariance matrix ($grid_size × $grid_size)."
    cached = T(bcf, t_end, grid_size)
    save_object(path, cached)
    logging && @info "Decomposition completed. Results saved to \"$path\"."

    return path
end


"""
    clear_covariance_cache()

Remove all cached `BCFEigen` objects.

If the cache directory does not exist, the function does nothing.
"""
function clear_covariance_cache()
    dir = "covariance_cache"
    isdir(dir) || return nothing

    # Remove all cached eigendecompositions.
    for entry in readdir(dir; join=true)
        rm(entry; recursive=true, force=true)
    end
    return nothing
end



function bcf_hash(bcf::BCF{N}) where {N}
    key = (
            N,
            bcf.Γ...,
            bcf.G...,
            (nameof(typeof(f.f)) for f in bcf.f_vector)...,
            (x for p in (f.params for f in bcf.f_vector) for x in p)...,
            (nameof(typeof(g.f)) for g in bcf.g_vector)...,
            (x for p in (g.params for g in bcf.g_vector) for x in p)...,
        )
    return hash(key)
end
