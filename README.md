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

```julia
using Pkg
Pkg.add(url="https://github.com/VladislavSukharnikov/MyLib.jl")
```

---

## Usage

First, load the project

```julia
using Humulus
```

A few examples of usage are given in folder examples.

First, one has to initialize all parameters required for the integration. One of the required objects is the bath-correlation function (BCF). While internal functions allow for definition of any type of BCF, the public constructors include squeezed-reservoirs, for instance, in example1.jl


```julia
bcf = let
    ω₀ = 5.0       # central (carrier) frequency
    Γ  = 1.0       # spectral half-width
    r  = 1.5       # squeezing parameter
    φ  = 0.0       # squeezing phase
    γ  = 1.0       # atom–field coupling strength
    
    one_mode_squeezed_bcf(ω₀, Γ, r, φ, γ)
end
```
this creates a BCF functor, that can be called as bcf(t,s) and has fields that store essential information about the BCF, for instance, bcf.Γ includes decay rates for all effective modes. The number of modes can be found as Humulus.n_modes(bcf).

The second object is the atomic parameters. The current code is implemented only for a two-level atom, and atomic parameters accept only the resonance frequency and initial (pure) state condition

```julia
atom_params = let
    # central (carrier) frequency
    ν₀ = 5.0

    # Initial atomic state:
    # |ψ⟩ = c_g |g⟩ + c_e |e⟩
    c_e = inv(√2)
    c_g = inv(√2) * exp(-1im * π/4)

    AtomParams(ν₀, c_g, c_e)
end
```
The output is the immutable struct.

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