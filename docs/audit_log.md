# Audit log — Bombadillo

When and why each decision was taken. Appended, never rewritten. Not read by default.

One entry per major decision or discovery:

```
## YYYY-MM-DD — <hash> — <title>
**Discovered / Decided:** what.
**Why:** the reasoning, including what was rejected.
**Cost:** what it broke, or left open.
```

**The hash names the commit that holds the work**, so the entry is written after it, and says
`uncommitted` until then. `dev-skills/logging` writes entries and fills in hashes. Entries above
the `<!-- hashes from here -->` line predate hashes and carry none.

---

## 2026-11-05 — `AGENT.md` created

- **What:** `AGENT.md` added, adapted from Hermetica's `AGENT.md`
  (github.com/patrickCNMartin/Hermetica). Behaviour kept: how to keep the file, the
  three-file split, principles 1–9, ponytail, auditing, "as diff".
- **Dropped as Hermetica-specific:** architecture invariants (hashing, artefacts,
  lifecycle, lock files, rendering, pipelines, storage), protocols.io rules, Python
  entry points, ruff, CI, containers and images, committed test fixtures, `.gitignore`
  rules, `detect-secrets`, the `~/.claude/settings.json` memory settings, `docs/status.md`.
- **Dropped:** everything about Nix. It is not available in this workspace; the "run" rule
  is the project-local Julia via `/project/.julia/env.sh` instead.
- **Adapted:** Python idioms rewritten for Julia (inner constructors, concrete field types,
  `Exception` subtypes with `showerror`, `ENV` instead of `os.getenv`); principle 3 now
  requires an explicit RNG/seed instead of `.env` config; the "as diff" tool is
  `git diff --no-index` because VS Code is not installed here.
- **Added for this repo:** the memory-cap rules, the core/legacy split ("legacy is not
  dead while the port is open"), and the `test/` vs `tests/` distinction.
- **Open:** no `.gitignore` exists, so `AGENT.md`, `docs/`, `src/core/` and `test/` are
  untracked and will be added by default on the next commit.

## 2026-11-05 — `.gitignore` added, README rewritten

- **`.gitignore`:** `.DS_Store` and `LocalPreferences.toml`. The two tracked `.DS_Store`
  files were untracked with `git rm --cached` (still on disk).
- **Why `LocalPreferences.toml` is ignored:** its settings are tuned to this workspace's
  memory cap, not needed on a normal machine. It stays on disk here and is still
  load-bearing locally.
- **Kept tracked on purpose:** `Manifest.toml` (root, already tracked), `test/Manifest.toml`,
  `test.gif`.
- **README:** corrected from "python module" to Julia; added install, plotting and test
  lines.

<!-- hashes from here -->

---

## 2026-11-12 — uncommitted — decisions are logged by `dev-skills/logging`, with their commit's hash

**Decided:** `dev-skills/logging` (tracked) writes audit-log entries and fills in their hashes.
The work is committed first, locally, staging only the paths the decision touched; the entry
follows in a commit of its own. `dev-skills/logging/pending.sh` lists entries whose hash is
`uncommitted`, missing, or no commit in this repo. AGENT.md runs it at the start of every session,
and the agent raises what it lists with the user before resolving anything. Nothing is pushed.

**Why:** a decision is taken before its commit exists, so the hash was left for later and entries carried none. Asked for.

**Cost:** a logged decision is two commits, the work and then the entry. Entries above `<!-- hashes from here -->` keep no hash; the check skips them.
