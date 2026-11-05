# Findings — Bombadillo

What was discovered: measurements, traps, probe results, superseded beliefs. Read when a
rule in `AGENT.md` looks arbitrary, or before re-investigating something. Each entry says
where the evidence lives.

---

## Workspace memory cap

- **The workspace has a ~2 GiB memory cap shared by every process**, including other agent
  sessions. Julia's defaults (one thread per visible core, parallel precompile, automatic
  precompile on `Pkg` operations) exceed it. Source: `/project/.julia/env.sh`.
- **`env.sh` sets a `ulimit -v` ceiling (default 1 300 000 KB)** so Julia aborts on its own
  instead of taking the workspace down. Override per run with `JULIA_VMEM_KB=…` — for
  example, a JET/AllocCheck run that needs more headroom.
- **`MALLOC_ARENA_MAX=2` is needed for that ceiling to work**: glibc reserves 64 MB of
  *virtual* memory per malloc arena, which uses up a `ulimit -v` budget on space that is
  never touched.
- **Extra threads cost memory, not just CPU.** 20 cores are visible; compute, GC and
  precompile image threads are all pinned to 1.
- **`JULIA_PKG_UNPACK_REGISTRY=true`** reads the General registry from unpacked files
  instead of holding the tarball in memory.
- **A precompile that dies with no Julia error is the ceiling.** Check this section before
  debugging it as a code error.

## Precompile workloads

- **DataFrames (root env) and JET, AllocCheck, GPUCompiler, LLVM (test env) run precompile
  workloads larger than the cap.** `precompile_workload = false` in the two
  `LocalPreferences.toml` files switches them off. The packages still work; their first
  calls are JIT-compiled and slower.
- **Preferences only apply to packages named in the active project.** GPUCompiler and LLVM
  are indirect dependencies of AllocCheck/JET, so `test/Project.toml` lists them under
  `[extras]` only so that their preferences take effect. Do not remove those entries as unused.

## Repo state

- **`LocalPreferences.toml` is gitignored**, so a fresh clone in this workspace lacks it.
  Recreate both files (root and `test/`) before precompiling, or the workloads above hit
  the cap.
- **`flake.nix` / `flake.lock` target machines with Nix.** Nix is not available in this
  workspace; the toolchain here is `/project/.julia` via `env.sh`.
- **The legacy repressilator in `src/cycle.jl` never oscillated.** Commits `5fbd739` and
  `07b2821` ("differential decay: doesn't work either — the cycling is just not happening")
  record the failed attempts. This is the motivation for the array-based core in
  `src/core/`.
- **`tmp/scan_rep.jl`** is a parameter scan (K, β, Hill n, k_off, 3 seeds) of the core
  repressilator that measures oscillation by autocorrelation (first negative minimum, then
  the next maximum as the period). Its results are not recorded here; rerun it before
  quoting them.
