# Humulus.jl

**Humulus** is a Julia implementation of the **Hierarchy of Pure States (HOPS)** method for simulating non-Markovian open quantum systems interacting with nonstationary Gaussian reservoirs. The name comes from *Humulus lupulus*, Latin for “hops.”

The repository currently contains implementations of **HOPS** and the associated **Hierarchy of Master Equations (HME)** for a two-level atom coupled to a squeezed reservoir.

This implementation accompanies the paper:

> V. Sukharnikov, S. Chuchurka, and F. Schlawin,
> *Non-Markovian dynamics in Nonstationary Gaussian Baths: A Hierarchy of Pure States Approach*,
> Physical Review Research (2026).
> DOI: https://doi.org/10.1103/yt37-s9hz

---

## Installation

For now, the code has not yet been unified into a single module. To set it up, include the following file:

```julia
include("src/Humulus.jl")
```

This file imports the required packages, includes all necessary files, and defines constants.

---

## Source Files

| File | Description |
|------|-------------|
| `input.jl` | Definition of model and simulation parameters. |
| `bcf.jl` | Construction of bath correlation functions. |
| `fock_space.jl` | Construction of the truncated pseudo-Fock space.  |
| `setup.jl` | Project setup. |
| `solver_params.jl` | Assembly of parameters required by the numerical solvers. |
| `hme/hme.jl` | Solver for the Hierarchy of Master Equations. |
| `hops/hops.jl` | Solver for the Hierarchy of Pure States. |
| `hops/noise.jl` | Generation of stochastic noise. |

---

## Examples

First, load the project

```julia
include("src/setup.jl")
```

The example of usage is provided in

```julia
include("examples/example1.jl")
```

<!-- ---


## Documentation

A detailed description is provided in **documentation.pdf**. -->

---

## Citation

If you use this software in academic work, please cite

```bibtex
@article{sukharnikov2026,
  author = {Sukharnikov, Vladislav and Chuchurka, Stasis and Schlawin, Frank},
  title = {Non-Markovian dynamics in Nonstationary Gaussian Baths: A Hierarchy of Pure States Approach},
  journal = {Physical Review Research},
  year = {2026},
  doi = {10.1103/yt37-s9hz}
}
```