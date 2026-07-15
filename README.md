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
First, load the package:

```julia
using Humulus
```

Several usage examples are provided in the `examples` directory.

Before running a simulation, all parameters required for the integration must be initialized. One of the required objects is the bath-correlation function (BCF). Although the internal API supports arbitrary BCF definitions, the public constructors currently include squeezed reservoirs. For example, in `example1.jl`:

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

This creates a BCF functor that can be called as `bcf(t, s)`. The returned object also stores parameters associated with the BCF. For example, `bcf.Γ` contains the decay rates of all effective modes. The number of modes is given by `Humulus.n_modes(bcf)`. Another available public constructor is `three_mode_squeezed_bcf`.

The second required object specifies the atomic parameters. The current implementation supports only a two-level atom. The atomic parameters consist of the resonance frequency and the initial pure state:

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

This returns an immutable `AtomParams` object.

The third required object defines the time grid used for the numerical integration:

```julia
grid_params = let
    t_end    = 20.0  # final simulation time
    n_save   = 250   # number of intervals between saved time points
    substeps = 4     # maximum internal integration substeps per save interval

    @info "Expected minimal number of steps is $(n_save * substeps)."

    GridParams(t_end, n_save, substeps)
end
```

The differential equations are solved using the `DifferentialEquations.jl` library, specifically the adaptive `Tsit5` solver. The `substeps` parameter determines the maximum allowed time step, `dt_max`, which is stored in `grid_params.dt_max` and computed as the spacing between consecutive saved time points divided by `substeps`.

Finally, the pseudo-Fock space truncation must be specified:

```julia
# Truncation of the pseudo-Fock basis
max_occupancy = 50
```

If the number of modes is $N$, the truncation can be specified either as an $N$-element tuple, an $N$-element vector, or a single integer. In the latter case, the basis includes all occupation-number configurations satisfying

$$
\sum_i n_i \le \text{max_occupancy}.
$$

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