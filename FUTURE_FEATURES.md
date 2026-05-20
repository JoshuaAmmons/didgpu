# Future features — companion tools for the DiD robustness workflow

These candidates come from auditing a real applied-econometrics
project (a DAO governance event-study with a 74 MB panel) for the
methods researchers actually run *alongside* the modern DiD
estimators. didgpu already covers the hard, modern estimators — the
de Chaisemartin/D'Haultfoeuille, Callaway-Sant'Anna, fect, TestMechs,
and HonestDiD families. The gaps below are the *companion* tools that
currently force a researcher to stitch together four or five packages
in one script.

## What real usage already maps cleanly onto didgpu

| Method in the wild | Reference package | didgpu equivalent |
|---|---|---|
| de Chaisemartin DiD (main, subgroups, LOO) | `DIDmultiplegtDYN` | `didgpu()` (60x faster, bit-identical) |
| Callaway-Sant'Anna | `did` | `didgpu_cs()` |
| HonestDiD (RM + smoothness) | `HonestDiD` | `didgpu_honest_did()` |
| Joint pre-trend Wald | (built into the DiD pkgs) | `p_jointplacebo` |
| Leave-one-out | hand-rolled loops | `didgpu_loo()` |
| Subgroup analysis | hand-rolled loops | `didgpu_by()` |

## Gaps to fold in (priority order)

### 1. `didgpu_twfe()` — naive TWFE / OLS event-study baseline  (task #98, HIGH)
Two-way fixed-effects event study, distributed-lag form for
non-absorbing treatment:

    Y_it = alpha_i + lambda_t
           + sum_{k>=0} beta_k  D_{i, t-k}   (Effect_k)
           + sum_{j>=1} gamma_j D_{i, t+j}   (Placebo_j)
           + eps_it

with cluster-robust SEs. This is the universal "naive" comparison
reported in nearly every applied DiD paper, and the Liu-Makridis
model in the audited project (built ad hoc with
`fixest` + `sandwich` + `lmtest`). Output deliberately mirrors
`didgpu()`'s Effect_k / Placebo_k matrices so it drops into the same
table/plot/robustness pipeline. Validate point estimates bit-exact
vs `lm()`/`fixest`; cluster SEs vs a hand-computed CR1 reference.

### 2. Pre-trends TOST equivalence test  (task #99, HIGH)
Two-one-sided-tests equivalence on the pre-treatment placebo
coefficients: is each (or the joint) pre-period coefficient inside
[-delta, +delta]? The modern, better-powered complement to the joint
Wald (which can only *fail to reject* flat pre-trends). Extends the
existing Placebo_k machinery; increasingly expected by referees.

### 3. Windowed / subset joint placebo test  (task #100, MEDIUM, tiny)
Let the joint pre-trend Wald run on a chosen window of placebos
(e.g. Placebo_1..3, the "short-window" test) rather than all of them.
A trivial extension of machinery didgpu already has.

### 4. `didgpu_bacon()` — Goodman-Bacon decomposition  (task #101, MEDIUM)
Decompose the TWFE estimate into its 2x2 sub-comparisons and weights,
for *absorbing* staggered designs. Not applicable to non-absorbing
treatments (like the audited project), but a standard "reviewer asks
for it" diagnostic for the large class of users with absorbing
treatment. The de Chaisemartin/CS framework is the modern
replacement, so this is lower priority — included for completeness.

## Why this matters

Folding these in lets a researcher run their *entire* DiD robustness
suite — modern estimator + naive baseline + pre-trend tests +
sensitivity + leave-one-out — through one fast, GPU-aware package,
instead of juggling DIDmultiplegtDYN + did + HonestDiD + fixest +
bacondecomp with hand-rolled glue. The audited project's
leave-one-out alone drops from ~1.5 days to minutes just by switching
its inner refit from `did_multiplegt_dyn` to `didgpu()`.
