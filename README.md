# FESOM_F

FESOM_F is a clean-architecture **Fortran** reimplementation of the FESOM2 ocean + sea-ice model.
Its numerical kernels are transcribed line-for-line from FESOM2, so the model reproduces a reference
FESOM2 configuration **bit-for-bit**; the surrounding "plumbing" (data types, the parallel layer,
configuration, I/O) is redesigned for clarity and future flexibility. The guiding rule is
*"first reproduce exactly, then optimize."*

**Status:** v1 is feature-complete and byte-exact against FESOM2. Milestones M0–M10 are done and
tagged (`m0`, `m1`, `m2-mvp`, `m3` … `m10`): dynamical core, sea-ice EVP/mEVP, GM/Redi, KPP, zstar
ALE, TKE, multi-year production runs, Zarr output, and restart/checkpoint. The next milestone is open.

See **[IMPLEMENTATION.md](IMPLEMENTATION.md)** for the architecture and design decisions, and
**[TESTING.md](TESTING.md)** for how it is validated against FESOM2 and how it scales.

> **A note on the name.** The project is called **FESOM_F** in this documentation, but it lives in a
> repository still named `fesom3`, and you will see that working name in real identifiers: the
> environment-variable prefix `FESOM3_`, build directories like `build_intel_dp`, driver binaries
> `fesom_*`, and git tags `m0…m10`. Those are the actual names — use them verbatim. Only the prose
> uses "FESOM_F".

---

## 1. Requirements

FESOM_F currently builds and runs on the **Levante** supercomputer (DKRZ). Everything is provided
through environment modules, loaded automatically by `env.sh`:

| | Anchor (bit-identity) | Portability |
|---|---|---|
| Compiler | `intel-oneapi-compilers/2022.0.1` | `gcc/11.2.0` |
| MPI | `openmpi/4.1.2-intel-2021.5.0` | `openmpi/4.1.2-gcc-11.2.0` |
| netCDF | `netcdf-c/4.8.1` + `netcdf-fortran/4.5.3` | same |

- **netCDF-Fortran is required** (forcing and initial-condition reads); the build aborts if
  `nf-config` is not on the `PATH`.
- **liblz4 is optional** — it enables the LZ4 output compressor; without it the writer is
  uncompressed-only.
- The **Intel + double-precision** build is the *anchor*: its compiler flags are copied verbatim
  from FESOM2 v2.7.3, which is what makes bit-identical results possible. The GNU build is for
  portability testing.

---

## 2. Build

`configure.sh` sources `env.sh` for you (loading the modules and the Levante MPI workaround), runs
CMake, and builds:

```bash
./configure.sh --compiler intel --precision dp --clean --build   # anchor build
```

Binaries land in `build_intel_dp/bin/`, Fortran modules in `build_intel_dp/module/`. The portability
build is the same with `--compiler gnu` (→ `build_gnu_dp/`):

```bash
./configure.sh --compiler gnu --precision dp --clean --build
```

`configure.sh` options: `--compiler intel|gnu`, `--precision dp|sp|hp`, `--debug` (builds into a
separate `build_<compiler>_<precision>_debug/` so it never overwrites the anchor), `--clean`,
`--build`, `--jobs N`. The build directory name encodes the compiler and precision, so several builds
coexist.

> Build and run on a Levante node — `env.sh` only knows how to set up Levante and will stop on any
> other host.

---

## 3. Quick check

Run the self-test suite on a login node:

```bash
cd build_intel_dp && ctest --output-on-failure
```

For a stronger check that FESOM_F still reproduces FESOM2 bit-for-bit, run one byte-gate — it drives
the real FESOM2 timestep and the FESOM_F timestep on identical state and compares every value:

```bash
bash tools/run_step_gate.sh        # expect: max|Δ| = 0
```

What these checks mean — and the full validation story — is in [TESTING.md](TESTING.md).

---

## 4. Running a simulation

The production driver is **`fesom_lifecycle_native_mr`** (multi-rank, computes its own air–sea
fluxes from atmospheric forcing). A run is configured through `FESOM3_*` environment variables; the
example below mirrors `tools/run_lifecycle_2yr_freerun_f3_dist864.sbatch`:

```bash
source env.sh intel                         # modules + Levante MPI workaround

POOL=/pool/data/AWICM/FESOM2
export FESOM3_MESH_DIR="$POOL/MESHES_FESOM2.1/core2"
export FESOM3_IC_FILE="$POOL/INITIAL/phc3.0/phc3.0_winter.nc"
export FESOM3_FORCING_DIR="$POOL/FORCING/JRA55-do-v1.4.0"
export FESOM3_FORCING=JRA55
export FESOM3_RUNOFF_FILE="$POOL/FORCING/JRA55-do-v1.4.0/CORE2_runoff.nc"
export FESOM3_SSS_FILE="$POOL/FORCING/JRA55-do-v1.4.0/PHC2_salx.nc"
export FESOM3_START_CLOCK="0 1 1958"        # "timeofday day year"
export FESOM3_NSTEPS=35040                   # 2 model years at 48 steps/day

# physics (the default FESOM2 production configuration)
export FESOM3_WHICH_ALE=zstar
export FESOM3_FER_GM=1 FESOM3_REDI=1 FESOM3_MIX_TKE=1 FESOM3_SW_PENE=1 FESOM3_CHL_SWEENEY=1

srun -l -n 864 ./build_intel_dp/bin/fesom_lifecycle_native_mr > run.log 2>&1
```

**Where the knobs are defined.** The authoritative, always-current list of variables is the
environment block at the top of **`src/drivers/fesom_lifecycle_native_mr.F90`** (each is read with
`get_environment_variable`). The `tools/*.sbatch` and `tools/run_lifecycle_*.sh` scripts are working
examples. The most common variables:

| Variable | Meaning |
|---|---|
| `FESOM3_MESH_DIR` | mesh directory (must contain `dist_<NP>/`) |
| `FESOM3_IC_FILE` | initial T/S netCDF |
| `FESOM3_FORCING_DIR`, `FESOM3_FORCING` | atmospheric forcing directory and set (`JRA55`/`CORE2`) |
| `FESOM3_RUNOFF_FILE`, `FESOM3_SSS_FILE` | runoff and sea-surface-salinity restoring climatologies |
| `FESOM3_START_CLOCK` | cold-start clock, `"timeofday day year"` |
| `FESOM3_NSTEPS` | number of time steps (or use `FESOM3_RUN_LENGTH` + `FESOM3_RUN_UNIT`) |
| `FESOM3_WHICHEVP` | sea-ice solver: `0` = EVP, `1` = mEVP |
| `FESOM3_WHICH_ALE` | vertical coordinate: `linfs`, `zlevel`, `zstar` |
| `FESOM3_FER_GM`, `FESOM3_REDI` | Gent–McWilliams bolus advection, Redi isopycnal diffusion |
| `FESOM3_MIX_KPP`, `FESOM3_MIX_TKE`, `FESOM3_SW_PENE` | vertical mixing and shortwave penetration |
| `FESOM3_OUTPUT`, `FESOM3_OUTPUT_EVERY` | output directory (enables output) and step cadence |
| `FESOM3_RESTART`, `FESOM3_RESTART_LENGTH`, `FESOM3_RESTART_UNIT` | checkpoint directory and cadence |

> **Configuration today vs. intent.** The production drivers are configured by these `FESOM3_*`
> variables plus in-code defaults. A write-once configuration *module* (`mod_config`, reading a
> single `namelist.config`) already exists, but the lifecycle drivers do not read it yet;
> consolidating run configuration into one namelist file is the intended future interface.

**Restart.** Set `FESOM3_RESTART=<dir>` to write checkpoints (cadence via `FESOM3_RESTART_LENGTH` /
`FESOM3_RESTART_UNIT`). On the next launch, if a complete checkpoint is found under that directory
(a `restart.latest` pointer to a folder containing `checkpoint.json`), the run resumes from it
byte-exactly; otherwise it cold-starts. See [IMPLEMENTATION.md](IMPLEMENTATION.md) §I/O.

---

## 5. Using different meshes

A mesh is a directory holding the global mesh files plus one sub-directory of partition files per MPI
rank count:

```
core2/
  nod2d.out  elem2d.out  aux3d.out  nlvls.out  elvls.out  edges.out  edge_tri.out  edgenum.out
  dist_1/  dist_2/  dist_8/  dist_864/   ...      # one per rank count (partition files)
```

Point `FESOM3_MESH_DIR` at this directory. The driver **auto-selects `dist_<NP>/`** from the number
of MPI ranks (`srun -n <NP>`); a single-rank run reads the global `nod2d.out`/`elem2d.out` directly.
So switching meshes — or switching rank counts on the same mesh — is just changing `FESOM3_MESH_DIR`
and `-n`, provided the matching `dist_<NP>/` exists.

Two meshes are used routinely:

- **pi** — a small test mesh (~3k nodes), shipped with the FESOM2 oracle at
  `/home/a/a270088/port2/fesom2/tests/data/MESHES/pi` (with `dist_1/`, `dist_2/`, `dist_8/`). This is
  the default when `FESOM3_MESH_DIR` is unset, and what the byte-gates and `ctest` use — it runs in
  seconds on a login node.
- **CORE2** — the global ~1° production mesh (242 917 nodes) at
  `/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2`, with partitions up to `dist_864`.

To run on a new mesh, create its `dist_<NP>/` partition for your target rank count (FESOM2's mesh
partitioner produces these), then point `FESOM3_MESH_DIR` at the mesh directory.

---

## 6. Output and reading results

Set `FESOM3_OUTPUT=<dir>` to enable output. Variables and writer options are read from `namelist.io`
in the run directory (see `config/namelist.io`); without it, a default variable set is written at the
`FESOM3_OUTPUT_EVERY` step cadence. Output is written as **Zarr** stores — one per variable per
period:

```
<dir>/ssh.fesom.1958.zarr
<dir>/temp.fesom.1958.zarr
<dir>/fesom.mesh.diag.zarr        # mesh geometry / connectivity (UGRID-style)
```

The stores are partition-independent (the same regardless of rank count) and open directly with
xarray:

```python
import xarray as xr
ds = xr.open_zarr("ssh.fesom.1958.zarr")   # CF time axis, embedded lon/lat
```

---

## 7. Repository layout

```
src/            Fortran source
  params/         precision, constants, configuration            (6 modules)
  types/          mesh / partition / dynamics / tracers / ice    (5)
  infra/          MPI halo layer, partitioning, clock, timers    (9)
  mesh/           mesh read, geometry, rotation                  (4)
  oce/            ocean physics kernels                          (24)
  ice/            sea-ice dynamics, thermodynamics, coupling     (6)
  forcing/        atmospheric forcing read + bulk fluxes         (3)
  io/             Zarr output, decomposition, restart            (8)
  step/           timestep orchestration                         (2)
  drivers/        program entry points + byte-gate drivers       (32)
cmake/          compiler-flag and helper modules
config/         namelist.io (output configuration)
env/            per-host module/setup files (levante.dkrz.de)
test/           ctest self-tests
tools/          byte-gate scripts, oracle runner, diff tools, sbatch jobs
docs/           HANDOFF, LESSONS, plans (internal working logs)
configure.sh    build helper          env.sh   environment loader
```

---

## 8. Documentation map

- **[IMPLEMENTATION.md](IMPLEMENTATION.md)** — architecture, design decisions, differences from FESOM2.
- **[TESTING.md](TESTING.md)** — validation approach, byte-gate methodology, results, scalability.
- **`docs/HANDOFF.md`** — the detailed internal development log (milestone-by-milestone state).
- **`docs/LESSONS.md`** — accumulated gotchas and subtle reproducibility traps.
- **`docs/plans/`** — per-milestone implementation plans (and `completed/`).
