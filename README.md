# Humulus.jl

**Humulus** is a Julia implementation of the **Hierarchy of Pure States (HOPS)** method for simulating non-Markovian open quantum systems interacting with nonstationary Gaussian reservoirs, characterized by the bath-correlation function of the form

$$
\alpha(t,s) = \sum_{j=1}^{N} G_j^2 e^{-\Gamma_j |t-s|} f_j(t) g_j^*(s).
$$

The name comes from *Humulus lupulus*, Latin for “hops.”

The repository currently contains implementations of **HOPS** and the associated **Hierarchy of Master Equations (HME)** for a two-level atom coupled to a squeezed reservoir.


This implementation accompanies the open-access paper:

> V. Sukharnikov, S. Chuchurka, and F. Schlawin,
> *Non-Markovian dynamics in Nonstationary Gaussian Baths: A Hierarchy of Pure States Approach*,
> Physical Review Research (2026).
> DOI: https://doi.org/10.1103/yt37-s9hz



---

## Installation

Install with Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/VladislavSukharnikov/Humulus")
```

Alternatively, install from a local clone by cloning the repository:

```bash
git clone https://github.com/VladislavSukharnikov/Humulus.git
cd Humulus
```

Then start Julia with the project environment and install the dependencies:

```bash
julia --project=.
```

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

---

## Usage
First, load the package:

```julia
using Humulus
```

Usage example is provided in the `examples` directory.

Before running a simulation, all parameters required for the integration must be initialized. One of the required objects is the bath-correlation function (BCF). Although the internal API supports arbitrary BCF definitions, the public constructors currently include squeezed reservoirs. For example, in `examples/example.jl`:

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

This creates a BCF functor that can be called as `bcf(t, s)`. The returned object also stores parameters associated with the BCF. For example, `bcf.Γ` contains the decay rates of all effective modes. Another available public constructor is `three_mode_squeezed_bcf`.

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
    n_save   = 500   # number of intervals between saved time points
    substeps = 5     # maximum internal integration substeps per save interval

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
\sum\_i n\_i \le \text{max}\textunderscore\text{occupancy}.
$$


Now the non-Markovian problem can be solved in one of two ways.

### Solving HME

To solve the hierarchy of master equations (HME), run

```julia
ρ_s = solve_hme(
    grid_params,
    bcf,
    atom_params,
    max_occupancy,
)
```

This returns the reduced system density matrix `ρ_s` with dimensions `ρ_s[i, j, t]`, where `i` and `j` label the atomic states (ground and excited), and `t` indexes the saved time points.

### Solving HOPS

Alternatively, the hierarchy of pure states (HOPS) can be solved by specifying the number of trajectories:

```julia
n_trajectories = 1
out = solve_hops(
    grid_params,
    bcf,
    atom_params,
    max_occupancy;
    n_trajectories=n_trajectories,
    clear_cache=false,
    show_progress=true,
    logging=true,
);
```

The HOPS implementation caches the Cholesky decomposition of the noise covariance matrix. The keyword argument `clear_cache` determines whether the cached data is deleted after the computation completes. By default, `clear_cache=true`. To manually clear the cache at any time, call

```julia

Humulus.clear_cache()
```

The cache is intended primarily for distributed computations. In such settings, sending large matrices to workers causes significant serialization overhead. Instead, the decomposition is stored in a `.jld2` file, which each worker can load independently. This approach is designed for computational clusters with a shared filesystem, such as the DESY Maxwell Cluster.

The returned object `out` is an `ArrayPartition` containing two arrays. The first stores the accumulated density matrix over all trajectories, while the second contains the total number of trajectories. The physical density matrix is therefore obtained as

```julia
ρ_s = out.x[1] ./ out.x[2]
```

## Parallel execution
The package supports parallel computation of stochastic trajectories using `Distributed` library. First, initialize the worker processes:

```julia
using Distributed

n_workers = Sys.CPU_THREADS
addprocs(n_workers)

@everywhere using Humulus
```

Then pass the workers to `solve_hops` via the `workers` keyword argument:

```julia
n_trajectories = 100

out = solve_hops(
    grid_params,
    bcf,
    atom_params,
    max_occupancy,
    n_trajectories;
    clear_cache = true,
    show_progress = false,
    workers = workers(),
)
```

In this example, each of the `n_workers` worker processes computes `n_trajectories` independent trajectories. Therefore, the expected total number of trajectories is

```julia
n_workers * n_trajectories
```

and `out.x[2]` should equal this value.

For long-running calculations, it is often preferable to split the computation into multiple batches. This makes the calculation more robust: if one or more workers disconnect during execution, the results from completed batches are preserved. The `@batched` macro provides this functionality:

```julia
n_trajectories = 100
n_batches = 100

out = @batched n_batches solve_hops(
    grid_params,
    bcf,
    atom_params,
    max_occupancy;
    n_trajectories = n_trajectories,
    clear_cache = false,
    show_progress = false,
    logging = false,
    workers = workers(),
)
```

In this case, each worker computes `n_trajectories` trajectories in each of the `n_batches` batches. Therefore, the expected total number of trajectories is

```julia
n_workers * n_trajectories * n_batches
```

and `out.x[2]` should again equal this value upon successful completion.

---

## Citation

If you use this software in academic work, please cite

```bibtex
@article{sukharnikov2026,
  author = {Sukharnikov, Vladislav and Chuchurka, Stasis and Schlawin, Frank},
  title = {{Non-Markovian dynamics in Nonstationary Gaussian Baths: A Hierarchy of Pure States Approach}},
  journal = {Physical Review Research},
  year = {2026},
  doi = {10.1103/yt37-s9hz}
}
```
