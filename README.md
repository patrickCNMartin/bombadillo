# Bombadillo

A Julia package for generating synthetic, procedurally built spatial omics data: a tissue
of cells, gene regulatory networks driving each cell's state, and counts sampled from that
state. Every dataset is reproducible from its parameters and seed.

Work in progress: an array-based simulation core (`src/core/`) is replacing the original
generators.

```julia
using Pkg; Pkg.develop(path = ".")
using Bombadillo
using Plots   # optional, enables the plotting functions
```

Tests: `julia --project=test test/runtests.jl`

## Contributing

**Audit log.** Decisions taken with a coding agent are logged by the agent with
[dev-skills/logging](dev-skills/logging/SKILL.md): one file per decision in `docs/decisions/`, committed
with its work. Decisions you take yourself go in the pull request description.

- Read the log: `sh dev-skills/logging/log.sh`
- Don't edit `docs/audit_log.md`, and don't rename or delete a record.
