using Humulus
using StaticArrays
using LinearAlgebra
using JLD2
using RecursiveArrayTools
using BenchmarkTools

using Humulus:
    BCF,
    random_bcf,
    n_modes,
    phasecomb,
    FuncWrapper,
    FockSpace,
    BCFEigen,
    sampler_from_cache,
    HME,
    HOPS, 
    create_solver_params