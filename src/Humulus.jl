module Humulus

    # =============================================================================
    # Dependencies. 
    # =============================================================================

    # Differential equation solver
    using OrdinaryDiffEq: ODEProblem, ODEFunction, solve, Tsit5

    # Data structures and utilities
    using StaticArrays: SVector, MVector, setindex
    using RecursiveArrayTools: ArrayPartition
    using LinearAlgebra: mul!, eigen!, Hermitian

    # Input/output/caching
    using JLD2: save_object, load_object

    # Parallel computation
    using Distributed: addprocs, rmprocs, pmap, WorkerPool


    # =============================================================================
    # Constants. 
    # =============================================================================

    const invroot2 = inv(sqrt(2.0))


    # =============================================================================
    # Source files. 
    # =============================================================================

    include("bcf/types.jl")
    include("bcf/functions.jl")
    include("bcf/models.jl")
    include("bcf/utils.jl")

    include("fock/types.jl")
    include("fock/constructors.jl")
    include("fock/basis.jl")

    
    include("params/time_grid.jl")
    include("params/atom_params.jl")
    include("params/grid_params.jl")
    include("params/solver_params.jl")

    include("hme/types.jl")
    include("hme/rhs.jl")
    include("hme/init_cond.jl")
    include("hme/solver.jl")


    include("noise/types.jl")
    include("noise/constructors.jl")
    include("noise/sampling.jl")
    include("noise/interpolation.jl")
    include("noise/cache_io.jl")

    include("hops/types.jl")
    include("hops/rhs.jl")
    include("hops/init_cond.jl")
    include("hops/solver.jl")

    include("bloch_vector.jl")

    # =============================================================================
    # Exports.
    # =============================================================================

    export  one_mode_squeezed_bcf, 
            three_mode_squeezed_bcf,
            AtomParams,
            GridParams,
            solve_hme,
            solve_hops,
            bloch_vector,
            @batched
    nothing
end 