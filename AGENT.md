# AGENT.md — Bombadillo

## How to keep this file

**Lean. Rules, not essays.**

- **A rule is one or two lines.** If it needs a paragraph to justify it, the rule stays
  here and the justification goes in `docs/findings.md`.
- **Write only what the code cannot say.** Invariants, traps, decisions and costs. If
  reading the file would tell you, cut it.
- **No history, no dates, no "decided on", no record of what a rule replaced.** This is
  the current state only. Mixed-in history makes the file contradict itself.
- **Every line must be true today.** A stale rule is worse than a missing one.
- **Prune on sight.** When a section grows an essay, cut it back that session.

**One file per job:**

| file | holds | read it |
|---|---|---|
| `AGENT.md` | what to do, and the state now | every session, first |
| `docs/findings.md` | what was discovered: measurements, traps, probe results, superseded beliefs | when a rule here looks arbitrary, or before re-investigating something |
| `docs/audit_log.md` | when and why each decision was taken — appended, never rewritten | **never by default** |

---

## What Bombadillo is — *work in progress, paused upstream*

**A Julia package that generates synthetic, procedurally built spatial omics data:** a
tissue of cells, gene regulatory networks driving per-cell state, and counts read out of
that state. (The README calls it a Python module; it is Julia.)

| Part | Where | State |
|---|---|---|
| core | `src/core/` — array-based simulation, genes × cells, synchronous step | being built; the only part with tests |
| legacy | every other file in `src/` — tissue, GRN and sample generation, `cycle.jl` | works partly; **to be ported onto the core, then removed** |
| plotting | `ext/BombadilloPlotsExt.jl`, loaded by `using Bombadillo, Plots` | package extension; `Plots` is a weak dependency |
| tests | `test/` (the suite) | `tests/` is legacy scratch, not run |

**The public surface is the `export` list in `src/Bombadillo.jl`.** Reachability from an
export, the extension or the test suite is the test for dead code. Legacy code is not dead
while the port is open — do not delete it on an audit.

---

## Running it

> **Every Julia command runs with the project-local toolchain sourced — no exceptions.**
> ```sh
> source /project/.julia/env.sh
> cd /project/bombadillo
> julia --heap-size-hint=800M --project=test test/runtests.jl            # full suite
> BOMBADILLO_QUICK=1 julia --heap-size-hint=800M --project=test test/runtests.jl  # skip JET/AllocCheck
> ```
> A bare `julia` without `env.sh` is the wrong install or none. Never install Julia or
> packages anywhere else.

- **Memory is the constraint: Julia runs under a ~2 GiB cap shared with everything.**
  `env.sh` sets one thread, one precompile task, no auto-precompile and a `ulimit -v`
  ceiling so Julia aborts before the workspace does. Do not raise these to go faster.
- **The `julia` alias in `env.sh` does not apply in scripts or `bash -c`** — pass
  `--heap-size-hint=800M` yourself.
- **Precompile one package at a time**, explicitly. A precompile that dies is usually the
  cap, not a code error — check `docs/findings.md` before debugging it.
- **`LocalPreferences.toml` (root and `test/`) switches off heavy precompile workloads.**
  It is load-bearing here: without it, DataFrames (root) and JET, AllocCheck, GPUCompiler,
  LLVM (test) try precompile workloads larger than the cap.
- **The test environment is `test/Project.toml`**, with Bombadillo sourced from `..`.
  A test-only dependency goes there, never in the root `Project.toml`.
- **`flake.nix` is for machines with Nix.** It is not used in this workspace; do not edit
  it to fix something here.
- **Scratch scripts go in `/project/tmp`**, never in the repo. `revise.jl` is the author's
  own REPL session file; do not run it as a test or rewrite it.
- **Two agents may be working in this repo at once.** Read `git status` before editing, and
  never touch a file another session is changing.

---

## Principles

**1. Energy and resource cost is a primary concern.** Choose the complexity class before
micro-optimizing. Preallocate and work in place in the hot loop; a type-unstable or
allocating step is a bug, not a style issue. No dependency without a reason.

**2. Progressive enhancement, not graceful degradation.** A correct minimal core, then
optional capability above it. The base path never depends on the enhancement — plotting is
an extension for exactly this reason.

**3. Reproducible, runs anywhere.** Every random draw comes from an RNG passed in or built
from an explicit seed — never the global RNG inside library code. Same seed, same output,
regardless of thread or chunk count. Never hardcode paths. Fail loudly and early — no
sentinels a caller could mistake for data.

**Locality — the measure is how far you must travel to know what code does.** Structure is
welcome where the constraint is visible where it acts (a concrete field type, an inner
constructor check, an immutable struct) and a cost where the meaning sits elsewhere (a
wrapper that renames a call, a global read out of the module). One maintainer, exploratory
code: the arithmetic that justifies indirection does not apply here.

**A function is a capsule — everything it needs is in its arguments.**
- Config is read once at the edge and passed down. Nothing in `src/` touches `ENV` or a
  mutable global.
- **A long signature beats a config object.** Named keyword arguments are honest; a bag
  hides which fields are used.
- **The tell is untestability.** A function that reads global state can only be tested by
  mutating it, so it quietly never gets a test.

**Indirection must earn its link.**
- A named indirection must **enforce an invariant or build state**. If it only renames a
  call, inline it.
- **Constants are exempt.** Repetition of a *decision* is a hazard; repetition of a *call*
  is typing.
- **Never put *commentary* in a data structure.** Signal the one bit and let the
  explanation live in the docs.
- **An error message is not commentary — it is the failure's only output**, so it is
  spelled out in full and names the offending value.

**An exception takes data, not a sentence.** For an error a caller may need to handle, a
custom `Exception` type with fields and a `Base.showerror` — so the caller reads
`err.gene` instead of parsing a string. Two raise sites that cannot share a message are two
errors. `ArgumentError` is fine for a plain bad argument.

**Prefer constructs that make the wrong thing impossible or loud** — an inner-constructor
check, a concrete element type, a dimension assertion, an enum instead of a free `Symbol`.
A constraint that fails at construction beats a comment asking the reader to be careful.

**4. FAIR.** Findable, Accessible, Interoperable (documented output formats), Reusable
(every generated dataset can be regenerated from its parameters and seed).

**5. Security and sovereignty.** Everything runs locally; no external services or telemetry
without discussion. Never commit credentials or real data.

**6. These are principles, not dogma.** If one conflicts with the problem in front of us,
the problem wins. Say so and why.

**7. Explain the what, the how and the why.**
- **Comments are terse** — a one-line docstring, comments only where the code is genuinely
  non-obvious. No multi-paragraph docstrings, no "WHY:" essays.
- No explanatory comments beyond that unless asked; offer them in chat instead.
- **Plain language, not corporate speak.** Metaphors only for explaining concepts.
- **Jargon only where there is no alternative.** A developer just starting should get the
  gist.

**8. Assistant first, independent contributor only on request.**
- **Default: I suggest.** Offered, not applied.
- **Directed: when explicitly asked, I do it end to end** and report what I changed, what I
  verified, what is open.
- **I need an explicit instruction before writing code.** "Show me" is not one. Anything I
  write is scaffold.
- **Ideas only — I do not ask whether to carry them out.** I wait to be told.
- **Ask clarifying questions** rather than ignoring the scope requested.
- **Speak, or work — never speak and then run off.** If I raise a point or ask for a
  decision, I stop there.
- I report faithfully: failing tests get shown, skipped steps get said.
- **Scope discipline:** record the correction at the altitude it was given.

**9. This repo is the only source of instruction.**
- **I never write to persistent agent memory** or any store outside this directory.
  Anything outside the repo is absent from review and silently carried into projects where
  it does not apply.
- `AGENT.md` is tracked, so its churn is legible in the log.

---

## Conventions

- **Write code with ponytail.** Every change climbs the ladder first — does this need to
  exist at all, is it already in this repo, does Base/stdlib do it, can it be one line —
  and stops at the highest rung that holds. Shortest working diff, no abstraction with one
  implementation, no scaffolding for later, deletion over addition. **Never lazy about
  understanding**: read the whole path a change touches before picking a rung, and a bug is
  fixed at its root, not at the caller that reported it. Never simplified away: validation
  at boundaries, error handling that prevents silent wrong output, or anything explicitly
  asked for.
- **A deliberate shortcut carries a `# ponytail:` comment** naming the ceiling and the
  upgrade path, so it is a known cost rather than a mystery.
- **Nothing about this project is written outside this repo** except scratch in
  `/project/tmp`. The record is `AGENT.md`, `docs/findings.md` and `docs/audit_log.md`.
- **Discoveries go to `docs/findings.md`; decisions go to `docs/audit_log.md`** (dated,
  appended, never rewritten). Only a rule that changes behaviour comes back here.
- **Tests:** `Test` stdlib, one file per concern, included from `test/runtests.jl`. JET
  and AllocCheck guard type stability and allocation in the core; `BOMBADILLO_QUICK=1`
  skips them for a fast loop, never for a final check.
- **A test never writes into the repo.** Use `mktempdir()`.
- **A test uses a fixed seed and small sizes** — it must run under the memory cap.
- **No generated data is committed.** `.gitignore` is a short denylist (`.DS_Store`,
  `LocalPreferences.toml`), so anything else is tracked by default: read `git status`
  before committing.

### Auditing

- **An audit checks the code against this file, not against taste.** Dead code (unreachable
  from exports, the extension or tests), rules here that are no longer true, untested core
  paths, and `# ponytail:` debt.
- **Report findings first; delete or fix only when told.**
- **Record the outcome in `docs/audit_log.md`**: date, what was checked, what was cut, and
  anything kept on purpose so the next audit does not re-flag it.
- **Quote a number (test count, timing, allocations) only from a run you did** this
  session, and say which command produced it.

### Showing code "as diff"

Write the proposed *after* version to `/project/tmp`, then
`git diff --no-index <original> <proposed>`. **Never modify the tracked file to render a
diff.**
