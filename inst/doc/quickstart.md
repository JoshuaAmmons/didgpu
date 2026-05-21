# didgpu quickstart

A 5-minute tour of using `didgpu` end-to-end. Assumes you've installed
the package (see `../../README.md`).

## 1. Simulate a panel

```r
library(didgpu)
p <- didgpu_simulate_panel(
  n_units      = 100L,
  n_periods    = 20L,
  frac_treated = 0.6,
  tau_profile  = c(0.5, 1.0, 1.2, 1.0, 0.8),  # true effect at k = 0, 1, ...
  sigma        = 0.4,
  seed         = 17L
)
head(p)
#>   unit period D          Y
#> 1    1      1 0  1.4015700
#> 2    1      2 0  1.7423457
#> 3    1      3 0  0.9819574
#> 4    1      4 0  1.5340345
#> 5    1      5 0  0.9974540
#> 6    1      6 1  1.5849837
```

The simulator stores the true DGP in `attr(p, "truth")` so you can
recover the ground-truth `F_g`, `tau_profile`, etc.

## 2. Verify didgpu agrees with the reference

If you have DIDmultiplegtDYN installed, `didgpu_compare()` runs both
backends and reports any disagreement:

```r
res <- didgpu_compare(p, "Y", "unit", "period", "D",
                      effects = 4L, placebo = 2L)
res$pass
#> [1] TRUE
res$report
#>     block      col max_abs_diff fails
#>   Effects Estimate     1.11e-16 FALSE
#>  Placebos Estimate     0.00e+00 FALSE
#>       ATE Estimate     1.11e-16 FALSE
```

Bit-identical estimates within machine epsilon.

## 3. Estimate with checkpointing

For a long-running bootstrap on real data, pass `checkpoint_dir`:

```r
fit <- didgpu(
  df = p,
  outcome   = "Y", group = "unit", time = "period", treatment = "D",
  effects        = 4L,
  placebo        = 2L,
  bootstrap_reps = 100L,
  seed           = 1L,
  checkpoint_dir = "checkpoints/run1",
  backend        = "r"
)
```

While it runs you'll see one log line per iter; each iter is saved to
`checkpoints/run1/cells/b####.rds` as soon as it completes. If the
process dies for any reason — OOM, power loss, you cancelled — just
re-invoke with the same call to resume:

```r
fit <- didgpu(p, "Y", "unit", "period", "D",
              effects = 4L, placebo = 2L, bootstrap_reps = 100L, seed = 1L,
              checkpoint_dir = "checkpoints/run1", backend = "r")
#> [didgpu] resuming checkpoints/run1: 47/101 cells already done
#> [didgpu] cell b=47    0.62s   (48/101 total)
#> ...
```

Resuming with the same seed and config produces identical aggregated
output to an uninterrupted run.

## 4. Inspect and plot

```r
print(fit)

# broom-style summaries
didgpu_tidy(fit)
#>        term   estimate  std.error statistic      p.value   conf.low conf.high    kind
#> 1  Effect_1 0.49364716 0.07217548  6.838716 8.001236e-12  0.3521858 0.6351086  effect
#> 2  Effect_2 0.91783542 0.08114627 11.310321 1.155555e-29  0.7587916 1.0768793  effect
#> ...

didgpu_glance(fit)
#>   n_effects n_placebos n_switchers_e1 n_obs_e1 p_jointeffects p_jointplacebo
#> 1         4          2             60      438   1.026517e-87   0.0156...
```

For a quick base-R event-study plot:

```r
plot(fit)   # dots = estimates, vertical bars = CI, vertical dashed line at horizon 0
```

Or use ggplot2 for richer styling:

```r
library(ggplot2)
ev <- didgpu_event_study_data(fit)
ggplot(ev, aes(event_time, estimate, colour = kind)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c(effect = "steelblue", placebo = "grey50")) +
  labs(x = "event time", y = "ATT", colour = NULL)
```

Standard accessors work like with `lm`:

```r
coef(fit)           # named numeric vector of estimates
confint(fit)        # CI matrix at the stored ci_level
vcov(fit)           # bootstrap covariance of (effects, placebos)
```

For a per-period (normalised) interpretation:

```r
didgpu(p, ..., normalized = TRUE)   # divides effect_k by k for binary
```

To allow group-specific linear trends (e.g. when units have heterogeneous
secular growth that the standard parallel-trends assumption can't handle):

```r
didgpu(p, ..., trends_lin = TRUE)   # FD the outcome, cumulative recovery
# ATE is suppressed under trends_lin (the identifying assumption doesn't
# pin down a single average); per-event-time effects are bit-identical
# to DIDmultiplegtDYN's trends_lin = TRUE output.
```

For sub-group analysis (`by =` in the reference):

```r
p$region <- ifelse(p$unit %% 2L == 0L, "north", "south")
didgpu_by(p, by_var = "region", outcome = "Y", group = "unit",
          time = "period", treatment = "D",
          effects = 4L, bootstrap_reps = 100L,
          checkpoint_dir = "checkpoints/by_region")
# each subgroup writes to a subdir; resumable per subgroup
```

To extend an existing checkpoint with more bootstrap reps:

```r
didgpu_bootstrap_more("checkpoints/run1", df = p, extra_reps = 100L)
```

## 5. Backend choice and speedup

```r
didgpu_backend_info()
#     backend available                                             notes
# 1 reference      TRUE                            DIDmultiplegtDYN 2.2.0
# 2         r      TRUE binary, no controls; effects + placebos supported
# 3       cpu     FALSE                             stub; needs Rcpp port
# 4      cuda     FALSE                   stub; needs CUDA toolkit + nvcc
```

For the binary case the `"r"` backend is ~14× to ~60× faster than the
`"reference"` backend (the gap widens with panel size). See
`../../README.md` for the benchmark table.

To switch back to the reference at any time:

```r
didgpu(p, ..., backend = "reference")
```

The result object has the same shape regardless of backend, so
downstream code keeps working.

## 6. Reading a finished checkpoint back

```r
chk <- didgpu_load_checkpoint("checkpoints/run1")
chk$meta            # the run's config
nrow(chk$manifest)  # how many cells were committed

# Reconstruct the aggregated fit from disk without rerunning.
cells  <- didgpu_aggregate_cells("checkpoints/run1")
length(cells)
```
