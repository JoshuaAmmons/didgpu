# A backend must never answer a question it cannot compute.
#
# Regression. .backend_cuda()'s compatibility guard checked only
# controls / weight / trends_nonparam, so every OTHER option passed
# straight through to a kernel that does not implement it and a wrong
# number came back SILENTLY -- no warning, no fallback.
#
# `normalized` was the damaging case. With a multivalued treatment,
# backend "cuda" returned the UNnormalised effects, which do not vary
# with dose, so a dose-response analysis looked exactly as though the
# treatment had been binarised. Against DIDmultiplegtDYN with
# normalized = TRUE on a 3-level dose:
#     r     |diff| 5.551e-17
#     cpu   |diff| 5.551e-17
#     cuda  |diff| 8.465e-01
# and backend = "auto" resolves to cuda whenever a GPU is present, so
# this was the DEFAULT path on a CUDA machine.
#
# The guard is now identical to .backend_cpu()'s. These tests compare
# each fast backend against the r-backend, which is the reference
# implementation for every one of these options.

.opt_panel <- function(seed = 11L, multivalued = TRUE) {
  set.seed(seed)
  nu <- 90L; np <- 12L
  ufe <- stats::rnorm(nu, 0, 1); tfe <- stats::rnorm(np, 0, 0.3)
  Fg <- rep(Inf, nu)
  tr <- sort(sample(seq_len(nu), 54L))
  Fg[tr] <- sample(4:9, 54L, replace = TRUE)
  dose <- rep(0, nu)
  dose[tr] <- if (multivalued) sample(c(1, 3, 5), length(tr), replace = TRUE) else 1
  g <- expand.grid(period = 1:np, unit = 1:nu)
  g <- g[order(g$unit, g$period), ]
  on <- g$period >= Fg[g$unit]
  g$D <- ifelse(on, dose[g$unit], 0)
  g$Y <- ufe[g$unit] + tfe[g$period] + 0.3 * g$D +
         stats::rnorm(nrow(g), 0, 0.4)
  g$w <- 1
  g[, c("unit", "period", "D", "Y", "w")]
}

.eff <- function(p, bk, ...) {
  f <- suppressMessages(suppressWarnings(
    didgpu(df = p, outcome = "Y", group = "unit", time = "period",
           treatment = "D", effects = 3L, placebo = 0L, bootstrap_reps = 0L,
           backend = bk, verbose = FALSE, ...)))
  as.numeric(f$results$Effects[, "Estimate"])
}

.backends <- function() {
  bk <- c("cpu")
  if (isTRUE(tryCatch(didgpu_has_cuda_support(), error = function(e) FALSE))) {
    bk <- c(bk, "cuda")
  }
  bk
}

test_that("normalized = TRUE on a multivalued treatment agrees across backends", {
  p <- .opt_panel(multivalued = TRUE)
  ref <- .eff(p, "r", normalized = TRUE)
  # Guard the test itself: normalisation must actually bite here, or the
  # comparison proves nothing.
  unnorm <- .eff(p, "r", normalized = FALSE)
  expect_gt(max(abs(ref - unnorm)), 1e-3)
  for (bk in .backends()) {
    expect_equal(.eff(p, bk, normalized = TRUE), ref, tolerance = 1e-12,
                 label = paste("backend", bk))
  }
})

test_that("every option in the compatibility guard falls back correctly", {
  p <- .opt_panel()
  opts <- list(
    normalized           = list(normalized = TRUE),
    trends_lin           = list(trends_lin = TRUE),
    same_switchers       = list(same_switchers = TRUE),
    only_never_switchers = list(only_never_switchers = TRUE),
    weight               = list(weight = "w")
  )
  for (nm in names(opts)) {
    ref <- tryCatch(do.call(.eff, c(list(p, "r"), opts[[nm]])),
                    error = function(e) NULL)
    if (is.null(ref)) next          # option not valid on this panel
    for (bk in .backends()) {
      got <- tryCatch(do.call(.eff, c(list(p, bk), opts[[nm]])),
                      error = function(e) NULL)
      expect_false(is.null(got), label = paste(nm, bk))
      if (!is.null(got)) {
        expect_equal(got, ref, tolerance = 1e-12,
                     label = paste(nm, "on backend", bk))
      }
    }
  }
})

test_that("dont_drop_larger_lower reaches the panel prep on every backend", {
  # CUDA called .prep_panel() without forwarding this flag, so it was
  # silently ignored there while being honoured on r and cpu.
  set.seed(4)
  p <- .opt_panel()
  # Make some units non-monotone so the flag actually changes the sample.
  u <- unique(p$unit)[1:20]
  p$D[p$unit %in% u & p$period >= 9] <- 0
  ref_on  <- .eff(p, "r", dont_drop_larger_lower = TRUE)
  ref_off <- .eff(p, "r", dont_drop_larger_lower = FALSE)
  # Only meaningful if the flag moves the r-backend answer.
  if (max(abs(ref_on - ref_off), na.rm = TRUE) > 1e-8) {
    for (bk in .backends()) {
      expect_equal(.eff(p, bk, dont_drop_larger_lower = TRUE), ref_on,
                   tolerance = 1e-12, label = paste("backend", bk))
    }
  } else {
    succeed()
  }
})

test_that("the treatment column is used, not just its support", {
  # Direct answer to "does didgpu binarise?". Under normalized = TRUE the
  # dose must change the estimate; if it did not, the values would be
  # getting collapsed to an indicator.
  p <- .opt_panel(multivalued = TRUE)
  pb <- p; pb$D <- as.numeric(pb$D > 0)
  expect_gt(max(abs(.eff(p, "r", normalized = TRUE) -
                    .eff(pb, "r", normalized = TRUE))), 1e-3)
})
