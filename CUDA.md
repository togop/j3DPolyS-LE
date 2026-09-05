# Optional CUDA Monte-Carlo backend

The Fortran engine can run the independent-trajectory Monte-Carlo loop on a GPU.
This is an **optional** build: the default `make build` / CMake configure stays
gfortran + MPI + HDF5 and does not require a CUDA toolkit.

## Build

```bash
# CPU only (default)
cmake -B cmake-build -S . -DUSE_CUDA=OFF
cmake --build cmake-build --target 3dpolys_le

# or
make build

# CUDA backend (needs nvcc)
cmake -B cmake-build-cuda -S . -DUSE_CUDA=ON
cmake --build cmake-build-cuda --target 3dpolys_le

# or
make build-cuda
```

`-DUSE_CUDA=ON` enables `enable_language(CUDA)`, compiles `sim_src/mc_cuda.cu`,
defines `USE_CUDA` for the Fortran ISO_C_BINDING wrapper, and links `cudart`.
If nvcc is not found the configure step fails rather than silently dropping CUDA.

## Runtime

| Flag | Effect |
|---|---|
| (default) | Use CUDA when the binary was built with `-DUSE_CUDA=ON` **and** a device is visible. |
| `--no-cuda` | Force the existing sequential `do_simulation` path. |
| `--cuda` | Request CUDA. If the binary was not built with `USE_CUDA`, this is a warning / no-op. |

If CUDA was compiled in but no device is present (or `--no-cuda`), the rank
falls back to sequential `do_simulation`. Each MPI rank still owns
`rank_Niter` trajectories. The GPU path launches **one thread per trajectory**
of that rank (batch size `T = rank_Niter`) and does **not** parallelize inside
a single polymer.

The log reports whether CUDA was compiled, whether a device was found, `T`,
and the method (`CUDA` vs `sequential`).

Init (helix / zigzag / `s=folder`) stays on the host. At each measurement the
device state is copied back and the existing `output()` writers produce
`config_*.out`, `dr_*.out`, `contact_*.out`, `Nlef_*.out` (and `process_*.out`
unbind events) so `analyse.f03` is unchanged. Snapshots are written in
**trajectory-major** order, same as the CPU loop. `hic3d` is not implemented
on the GPU.

## When it helps

One thread per independent replica. A single trajectory, or a tiny `Niter`,
will usually be **slower** than the CPU path (kernel launch + host I/O dominate).
The GPU path is aimed at many replicas per rank.

## RNG

Device RNG is **xoroshiro128+** (two `uint64` per trajectory), seeded from
`--seed` plus the trajectory index. Runs are reproducible on GPU but **not**
bit-identical to Fortran `random_number()`.

With `coherent_moves` (default **true** on CUDA) the *move type*
(tad / extrude / unbind / bind) is drawn from a launch-wide `hash_u01`
(splitmix64) counter. Per-move uniforms (monomer index, Metropolis accept,
Bernoulli `ikm`/`ikb`/`iku` draws) still come from xoroshiro. `pt = N / (3N + nfree[t])`
stays per trajectory because `nfree` differs. This is a valid Markov chain.

## Z-loop (+) swap fix

`trialmoveex` for a `+` leg that swaps with the occupant of `n+1` used
`config(2, n+1)` (bond `n+1 → n+2`). That is wrong: the bond being crossed is
`config(2, n)` (`n → n+1`). The same correction is applied on CPU and GPU.
The `−` branch was already correct.
