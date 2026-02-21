# GPU optimization design for `do_simulation`

This document describes how to accelerate the main model's `do_simulation` routine using GPUs (CUDA when available, or OpenACC), while **preserving the stochasticity** of the Monte Carlo trial moves.

## 1. Why replica parallelism (no change to move order)

The inner loop of `do_simulation` is:

```fortran
do v = 1, Ntrial
   r = randomnumber()
   if (r < pt)      call trialmovetad()
   elseif (r < 2*pt) call trialmoveex()
   ...
end do
```

Each trial move **reads and updates** shared state (`config`, `bittable`, `contact`, `dr`, `Nleffree`). The next move depends on the outcome of the previous one. So **within a single trajectory**, the trial sequence is inherently sequential and cannot be parallelized over `v` without changing the algorithm or introducing bias.

**Conclusion:** The only safe way to use the GPU without altering the Monte Carlo statistics is to run **multiple independent trajectories (replicas) in parallel**. Each replica has its own state and its own RNG stream; stochasticity is preserved because each trajectory remains a valid, independent realisation.

## 2. Algorithm: replica-parallel `do_simulation`

- **Replicas:** Run `N_replica` trajectories in parallel (e.g. one per GPU thread or block).
- **Per-replica state:** Each replica has its own copy of:
  - `config`, `bittable`, `contact`, `dr`, `boundary`, `loading_sites_factor`, `interaction_sites_state`
  - Scalar parameters: `L`, `Nchain`, `Nleffree`, `kb`, `ku`, `km`, `Ea`, `Ei`, etc.
- **Per-replica RNG:** Each replica must use an **independent** random stream. Use a dedicated RNG state per replica, seeded e.g. with `seed = base_seed + trajectory_id` (and optionally `+ rank * Niter` for MPI). Do **not** share one RNG across replicas.
- **I/O:** `output()` writes to host file units and uses global state. So:
  1. Run all replicas on the device for one “measurement interval” (all `Ninter` × `Ntrial` steps).
  2. After each measurement, leave the parallel region, copy the chosen replica states to the host (or keep them in host-accessible memory).
  3. On the host, for each replica in turn, set the model state and call `output()` (and optionally `log%info`), so that file layout and behaviour match the current one-trajectory-per-call design.

No change is made to the **order or logic** of trial moves inside a trajectory; only **which trajectories** are advanced at the same time changes (many in parallel instead of one after the other).

## 3. CUDA approach

- **Kernel design:** One trajectory per **thread block** (or per **thread** if state is small enough). Each block has its own copy of the model state in shared/global memory and its own RNG state (e.g. cuRAND).
- **RNG:** Use **cuRAND** (e.g. `curand_init(seed, 0, 0, &state)` per block/thread, then `curand_uniform_double()`). This gives one independent stream per replica and preserves stochasticity.
- **Code layout:** Either:
  - **CUDA Fortran** (nvfortran): Write a kernel that runs the same trial loop as the current Fortran `trialmovetad` / `trialmoveex` / `trialbound` / `trialunbound`, with device arrays and a device RNG, or
  - **C/C++ CUDA** with a small Fortran interface: Implement the trial loop in C/C++/CUDA and call it from Fortran; less natural for a Fortran-heavy codebase.
- **Data:** Flatten model state to contiguous arrays (e.g. `config(2, Nchain, N_replica)` in global memory). Copy in initial state per replica; after each measurement, copy out the slice for each replica and perform output on the host.
- **When “CUDA is available”:** Use runtime detection (e.g. `cudaGetDeviceCount` or a config flag) to choose the CUDA path; otherwise fall back to the current CPU loop.

## 4. OpenACC approach

- **Parallelism:** Use a single **parallel loop over replicas** (trajectories). Each gang/worker (or each vector element, depending on how you map “replica” to the OpenACC model) runs one full trajectory. The **inner** loop over `v` (trial moves) stays sequential within that replica.
- **Directives (conceptually):**
  - Copy in: replicate initial state for all replicas (e.g. `config`, `bittable`, etc.) to the device.
  - `!$acc parallel loop` over `replica = 1, N_replica`, with **private** model state per replica (or one large array with last index = replica).
  - Inside the loop: same `do_simulation` logic (burn-in, then measurement loop with `Ninter`×`Ntrial` trials and periodic output handling). No I/O inside the parallel region; only compute.
  - After each “measurement”, exit the parallel region (or use a serial region), copy back the state for all replicas, then on the host call `output()` (and logging) for each replica in order.
- **RNG:** Use a **per-replica** seed so that each OpenACC “thread” has its own stream. With **nvfortran** and `-acc`, `random_number()` in device code is implemented with one state per thread, so call `random_seed(put=seed_array)` at the start of each replica’s work (with `seed_array` derived from `trajectory_id`). With **gfortran** and OpenACC, device `random_number()` support may be limited; then use a small device-side RNG (e.g. a minimal LCG or xoroshiro) with a seed per replica.
- **Data:** Either:
  - **Array of derived types:** `type(PolymerModel), allocatable :: replicas(:)`. Some compilers (e.g. nvfortran) support deep copy of such arrays to the device; then each replica is one element. Or
  - **Flattened arrays:** e.g. `config_r(2, Nchain, N_replica)`, `bittable_r(14, bittable_t, N_replica)`, etc., and a version of the trial routines that take a replica index and operate on these arrays. This is more portable and often easier to tune.
- **Build:** Compile with OpenACC enabled (e.g. nvfortran `-acc`, or gfortran `-fopenacc`). Use a preprocessor or CMake option so that when OpenACC is disabled, the code falls back to the existing sequential trajectory loop.

## 5. OpenACC vs CUDA: performance and effort

| Aspect              | OpenACC                           | CUDA (hand-written)                    |
|---------------------|------------------------------------|----------------------------------------|
| **Parallelism**     | Replica parallelism only (same as above). | Same.                           |
| **Performance**     | Typically **~70–95%** of hand-tuned CUDA for this pattern (replica-parallel, little shared memory, modest per-thread work). | Maximum control and tuning.    |
| **Portability**     | One source for GPU (NVIDIA) and, in principle, other accelerators. | NVIDIA only (unless abstracted). |
| **Effort**          | Add directives and per-replica data layout; reuse existing trial routines where possible. | Rewrite or wrap trial loop in CUDA; more code. |
| **RNG**             | Compiler-provided or small device RNG per replica. | cuRAND per block/thread.       |

So: **OpenACC is a good first step** for GPU acceleration with less code change and acceptable performance; **CUDA is an option** when you need the last bit of performance or tighter integration with other CUDA libraries.

## 6. Preserving stochasticity – checklist

- [ ] **One RNG stream per trajectory:** Never share one RNG state across replicas.
- [ ] **Deterministic seed per trajectory:** e.g. `seed = base_seed + trajectory_id` (and same formula with/without GPU so that CPU and GPU runs are reproducible and comparable if desired).
- [ ] **No reordering of moves within a trajectory:** The sequence of trial moves (and the logic of `trialmovetad`, `trialmoveex`, `trialbound`, `trialunbound`) is unchanged; only multiple trajectories are run in parallel.
- [ ] **Same acceptance rule:** Metropolis (or current) acceptance with `exp(-dE)` and the same random number usage; no “batched” accept/reject that would change the statistics.

With this, the GPU run produces the **same statistical ensemble** as the CPU run (for the same seeds and parameters), up to floating-point order and compiler differences.

## 7. OpenACC vs CUDA: when to use which

- **Use OpenACC first** if you want one codebase, good GPU speedup (often 70–95% of hand-tuned CUDA for this replica-parallel pattern), and less implementation effort. Compile with `nvfortran -acc` (or `gfortran -fopenacc`); each replica gets its own RNG state via `random_seed` or compiler-provided per-thread state.
- **Use CUDA** when you need the last bit of performance, want to integrate with other CUDA libraries (e.g. cuRAND, cuFFT), or need fine control over blocks/threads and memory. Expect roughly 10–30% more performance than OpenACC for the same algorithm, at the cost of more and NVIDIA-specific code.

## 8. Suggested implementation order

1. **OpenACC replica path:** Add a “replica” data layout (flattened arrays or array of models) and a single `!$acc parallel loop` over trajectories; keep the existing trial routines, use per-replica RNG seed, and do output on the host after each measurement. This gives immediate GPU acceleration with minimal algorithm change.
2. **Optional CUDA path:** If needed, add a CUDA kernel (or C++/CUDA layer) that does the same replica-per-block (or per-thread) loop with cuRAND, and call it from Fortran when CUDA is available.
3. **Conditional build:** Use `ENABLE_OPENACC` (and optionally `ENABLE_CUDA`) in CMake so that the code compiles and runs correctly when GPU support is disabled.

## 9. Code changes added in this repo

- **`polymer_model.f03`**:
  - `advance_one_measurement_interval(self, Ninter)` runs one measurement interval (Ninter × Ntrial trial moves) with no I/O; restores/saves per-replica `rng_state` when allocated so replica-parallel runs keep independent RNG streams.
  - **RNG state:** Each model may hold allocatable `rng_state`; it is set in `init()` (capture current stream) and used in `advance_one_measurement_interval` to restore before and save after the interval.
  - **Phases:** `run_burnin_phase(self, burnin)`, `run_burnout_only_phase(self, burnout)`, `run_burnoutM_one_phase(self, Ninter)` run burn-in and burn-out steps with no I/O for use by the replica path.
- **`gpu_replica_mod.f03`**:
  - `run_measurement_interval_replicas(replicas, Ninter, use_gpu)`: one measurement interval for all replicas; with OpenACC and `use_gpu=.true.` the loop can run on the device.
  - `run_burnin_replicas`, `run_burnout_only_replicas`, `run_burnoutM_one_replicas`: same pattern for burn-in and burn-out phases.
  - `do_simulation_replicas(replicas, nrep, Ninter, Nmeas, burnin, burnout, burnoutM, use_gpu)`: full replica-parallel run (initial output, burn-in, main loop, burnout, burnoutM, erase); output is done on the host after each phase/measurement (measurement-major order in the output files).
- **`3dpolys_le.f03`**: When `use_gpu` is true and `rank_Niter > 0`, the program allocates `replicas(rank_Niter)`, sets a per-replica seed (from `random_seed_value` and `trajectory_i`), inits each replica, then calls `do_simulation_replicas(...)`. Otherwise it keeps the sequential one-trajectory-at-a-time loop.
- **`GPU_DESIGN.md`**: This design document.

## 10. Why you might see no GPU speedup (and what to do)

### Likely causes

1. **Offload not happening**  
   OpenACC compilers often **do not support** “deep copy” of **arrays of derived types with allocatable components**. So `copy(replicas)` or `enter data copyin(replicas)` may not actually move `PolymerModel` (with its allocatable `config`, `bittable`, etc.) to the device. The parallel loop may then run on the **host**, so you see no GPU speedup.  
   **Check:** Build with **nvfortran** `-acc -Minfo=accel` and look for “Generating copyin(replicas)” and “Offloading” / “GPU” messages. If you see “copy not supported” or no offload, the kernel is not running on the GPU.

2. **Data transfer dominating**  
   Even when offload works, **copying all replica state to the device and back** for every measurement interval can cost more than the compute. The code now uses a **data region** (`enter data copyin(replicas)` once, `update host(replicas)` only when writing output, `exit data copyout(replicas)` at the end) to reduce transfers. If your compiler does not support deep copy, this region may still not move data correctly.

3. **Small work per replica**  
   Each replica runs a **sequential** chain of trial moves. If `Ninter` × `Ntrial` is small or the number of replicas is small, the GPU may be underused and launch/transfer overhead can dominate.

4. **Compiler / runtime**  
   **gfortran** with `-fopenacc` has limited OpenACC support (e.g. derived types with allocatables, `random_number` on device). For GPU runs, **NVIDIA nvfortran** with `-acc` is much more reliable.

### What was changed to help

- **Data region:** Replica data is kept on the device across measurement intervals; `update host(replicas)` is used only before each output block, then `exit data copyout(replicas)` at the end. This reduces host–device copies compared to copying every kernel launch.
- **`present(replicas)`:** The parallel loops use `present(replicas)` so that, when the data region is active, the kernel uses device-resident data instead of triggering a copy each time.
- **`!$acc routine seq`:** The trial routines and the advance/burn-in/burn-out phase routines are marked so the compiler can generate **device code** for them (needed for nvfortran to run the loop on the GPU).

### Better way to utilize GPUs (reliable speedup)

For **predictable GPU utilization** and speedup, the robust approach is to **avoid derived types with allocatables** in device regions:

- **Flattened replica layout:**  
  Use contiguous arrays indexed by replica, e.g.  
  `config_r(2, Nchain, N_replica)`,  
  `bittable_r(14, bittable_t, N_replica)`,  
  `contact_r(3, Nchain, N_replica)`,  
  etc., and a **single device subroutine** that takes a replica index and performs one measurement interval (or one trial step) using these arrays. No deep copy is required; standard `copyin`/`copyout`/`present` work.

- **One data region for the run:**  
  `!$acc data copyin(config_r, bittable_r, ...) create(...)` at the start, run all measurement intervals (and burn-in/burn-out) on the device with `present(...)`, and `update host` only when writing output. That minimizes transfer and keeps compute on the GPU.

- **More replicas / longer runs:**  
  Run many replicas (e.g. tens to hundreds) so the GPU has enough parallel work, and use enough measurement steps so that compute clearly outweighs transfer and launch overhead.

Implementing the **flattened layout + device subroutine** is more invasive (duplicate or refactored trial logic that takes replica index and flat arrays) but is the path that works across compilers and gives real GPU speedup.

### Getting performance when GPU shows no speedup

- **Use CPU multi-core (OpenMP):** When the GPU path does not offload (or you run with `--no-gpu`), the replica loop is parallelized with **OpenMP** on the CPU. Set `OMP_NUM_THREADS` to the number of cores you want (e.g. `export OMP_NUM_THREADS=8`). Build with OpenMP available (CMake finds it automatically on most systems). You should see speedup when running many replicas (e.g. `Niter` ≥ number of threads).
- **Verify GPU offload:** Build with **nvfortran** and `-Minfo=accel` (e.g. set `CMAKE_Fortran_FLAGS` to include `-acc -Minfo=accel`). Inspect the compiler output for “Generating copyin(replicas)” and “Offloading” or “GPU”; if you see “copy not supported” or no offload messages, the kernel is running on the host.
- **Run more replicas:** With OpenMP, more replicas (e.g. 16–64) give better CPU utilization. With a working GPU path (e.g. after implementing the flattened layout), more replicas improve GPU occupancy.

## 11. Algorithm change: flattened replica layout (implemented)

To get **real GPU execution**, the code now uses a **flattened replica layout** when `use_gpu` is true:

- **Flat arrays:** State for all replicas is stored in plain arrays: `config_f(2, Nchain, nrep)`, `bittable_f(14, bt, nrep)`, `contact_f(3, Nchain, nrep)`, `dr_f(3, Nchain, nrep)`, `boundary_f`, `interaction_sites_state_f`, `loading_sites_factor_f`, `rng_state_f`, `Nleffree_f(nrep)`. No derived type with allocatables in device regions.
- **Module `gpu_flat_mod`:** Packs replicas into these arrays, runs `!$acc enter data copyin(...)`, then a single `!$acc parallel loop` over replica index calling `advance_one_measurement_interval_flat(irep)` (and similar for burn-in/burn-out). Device routines (`trialmovetad_flat`, `trialmoveex_flat`, `trialbound_flat`, `trialunbound_flat`, `unbound_all_flat`) operate on the `irep`-th slice of the flat arrays. After each phase/measurement, `update host` is used and state is unpacked for output.
- **When it runs:** The main program calls `do_simulation_replicas_flat` when `use_gpu` is true (and `rank_Niter > 0`). With OpenACC enabled (e.g. nvfortran `-acc`), the parallel loop and device routines run on the GPU.
- **Unbind log (unit 13):** On the device, `write(13,...)` is not performed; the unbind-event log is not written for the flat GPU path. Simulation results and other outputs are unchanged.

## 12. Would CUDA help?

- **OpenACC vs CUDA:** OpenACC with the **flattened layout** is sufficient for the compiler to offload the replica loop to the GPU. Once the kernel runs on the device (plain arrays + `routine seq`), you should see speedup. If you still see none, verify offload with nvfortran `-Minfo=accel`.
- **CUDA** (hand-written kernels in C/C++, called from Fortran via `iso_c_binding` or a C wrapper) can give **extra** performance (often 10–30%) by tuning block/grid sizes, using cuRAND, and avoiding any remaining compiler overhead. It does **not** fix the fundamental issue: the original “array of derived type” path does not offload; the **algorithmic** fix is the flat layout, which is now implemented with OpenACC.
- **Recommendation:** Build with **nvfortran** and **`-acc`** (OpenACC), run with `use_gpu` true so the flat path is used. If you then see GPU speedup, consider CUDA only if you need the last bit of performance.

## 13. Verifying GPU build and offload

### 1. Check that the build supports GPU (OpenACC)

- **At configure time:** Run CMake with OpenACC enabled and the NVIDIA HPC SDK Fortran compiler:
  ```bash
  export FC=nvfortran
  mkdir -p build && cd build
  cmake .. -DENABLE_OPENACC=ON
  ```
  Look for in the CMake output:
  - `OpenACC enabled: -acc -Minfo=accel` (or `-acc` if you removed the info flag)
  - `Fortran compiler: .../nvfortran` (or path to nvfortran)

- **If you use a different compiler:** Set `FC` to your Fortran compiler before running CMake. For **gfortran**, OpenACC is enabled with `-fopenacc` but GPU offload is often limited; for reliable GPU execution use **nvfortran**.

### 2. Verify offload with nvfortran and `-Minfo=accel`

The build is set up so that when `ENABLE_OPENACC=ON` and the compiler is NVHPC/PGI, the flags `-acc -Minfo=accel` are used. That makes the compiler print acceleration/offload information for every relevant construct.

- **Build and capture the report:**
  ```bash
  cd build
  make 2>&1 | tee build.log
  ```

- **In `build.log`, look for:**

  - **Offload is happening:** Lines like:
    - `gpu_flat_mod.f03: ... Accelerator kernel generated`
    - `Generating copyin(...)` or `Generating copy(...)` for the flat arrays (`config_f`, `bittable_f`, etc.)
    - `Generating present(...)` for the same arrays in the parallel loop
    - `Offloading` or `GPU` in the line describing the loop or kernel

  - **Routine (device) code:** Lines like:
    - `advance_one_measurement_interval_flat`, `trialmovetad_flat`, etc. with `Accelerator routine generated` or similar

  - **Offload is not happening:** Warnings or messages such as:
    - `copy not supported`, `variable in unstructured data reference`, or no “Accelerator kernel” / “Offloading” for the replica loop

- **Optional: only acceleration messages**  
  To reduce noise, you can build with only the accel info:
  ```bash
  make VERBOSE=1 2>&1 | grep -E 'accel|Accelerator|Offload|GPU|copyin|present'
  ```

### 3. Optional: force accel info without changing CMake

To try `-Minfo=accel` without editing CMake (e.g. if the default is only `-acc`):

```bash
export FC=nvfortran
cmake .. -DENABLE_OPENACC=ON -DCMAKE_Fortran_FLAGS="-acc -Minfo=accel"
make 2>&1 | tee build.log
```

Then inspect `build.log` as in step 2.

### 4. Runtime check that a GPU is used

- **Environment:** On a machine with an NVIDIA GPU, set:
  ```bash
  export ACC_DEVICE_TYPE=nvidia
  export ACC_DEVICE_NUM=0
  ```
  (Optional; nvfortran often defaults to the first NVIDIA GPU.)

- **Run a short test** (e.g. a small config with a few replicas). If the flat GPU path is used and offload works, you should see GPU utilization (e.g. with `nvidia-smi` in another terminal) and typically better performance than the CPU-only run with `--no-gpu`.
