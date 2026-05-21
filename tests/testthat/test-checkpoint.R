test_that("init_checkpoint creates expected files and refuses double-init", {
  cdir <- tempfile("didgpu_test_")
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)

  meta <- list(
    panel_hash = "deadbeef", seed = 1L,
    bootstrap_reps = 5L, effects = 2L, placebo = 1L,
    outcome = "Y", group = "unit", time = "period", treatment = "D",
    package_version = "0.0.0.dev"
  )
  expect_invisible(didgpu_init_checkpoint(cdir, meta))
  expect_true(file.exists(file.path(cdir, "meta.json")))
  expect_true(file.exists(file.path(cdir, "manifest.csv")))
  expect_true(dir.exists(file.path(cdir, "cells")))

  expect_error(didgpu_init_checkpoint(cdir, meta),
               "manifest.csv already exists")

  # force = TRUE wipes and re-inits.
  expect_invisible(didgpu_init_checkpoint(cdir, meta, force = TRUE))
})

test_that("manifest round-trips after writing cells", {
  cdir <- tempfile("didgpu_test_")
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)
  meta <- list(panel_hash = "h", seed = 1L, bootstrap_reps = 4L,
               effects = 2L, placebo = 0L, outcome = "Y", group = "g",
               time = "t", treatment = "D", package_version = "x")
  didgpu_init_checkpoint(cdir, meta)

  for (b in c(0L, 1L, 3L)) {
    didgpu:::.save_cell(cdir, b = b,
                        value = list(coef = c(a = b * 0.1)),
                        wall_seconds = 0.01)
  }
  chk <- didgpu_load_checkpoint(cdir)
  expect_equal(nrow(chk$manifest), 3L)
  expect_equal(sort(chk$manifest$b), c(0L, 1L, 3L))

  # todo set = expected reps minus done.
  expect_equal(didgpu:::.cells_todo(chk$manifest, 4L), c(2L, 4L))

  agg <- didgpu_aggregate_cells(cdir)
  expect_equal(sort(as.integer(names(agg))), c(0L, 1L, 3L))
  expect_equal(agg[["1"]]$coef, c(a = 0.1))
})

test_that("panel hash is order-independent in row order but row-content-sensitive", {
  p <- small_panel()
  h1 <- didgpu:::.panel_hash(p, "Y", "unit", "period", "D")
  h2 <- didgpu:::.panel_hash(p[sample(nrow(p)), ], "Y", "unit", "period", "D")
  expect_equal(h1, h2)
  p2 <- p; p2$Y[1L] <- p2$Y[1L] + 1
  h3 <- didgpu:::.panel_hash(p2, "Y", "unit", "period", "D")
  expect_false(identical(h1, h3))
})

test_that("didgpu_resume forwards all stored args (incl. normalized, switchers)", {
  # Start a run that uses several newer args. Stop after a couple cells
  # by truncating manifest. Resume via didgpu_resume(dir, df). The
  # aggregated result must match a fresh-from-scratch run with the same
  # args; that proves resume picked up every stored arg.
  p <- small_panel()
  cdir <- tempfile("didgpu_resume_")
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)

  # Run the full thing once to produce the "truth" aggregate.
  full <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 2L, placebo = 0L,
                  switchers = "in", normalized = TRUE,
                  bootstrap_reps = 4L, seed = 7L,
                  checkpoint_dir = cdir,
                  backend = "r", verbose = FALSE)

  # Now drop the last two committed cells and the manifest rows for
  # them, then ask didgpu_resume to top them back up.
  cdir2 <- tempfile("didgpu_resume2_")
  dir.create(cdir2)
  on.exit(unlink(cdir2, recursive = TRUE), add = TRUE)
  file.copy(file.path(cdir, "meta.json"),     file.path(cdir2, "meta.json"))
  file.copy(file.path(cdir, "manifest.csv"),  file.path(cdir2, "manifest.csv"))
  dir.create(file.path(cdir2, "cells"))
  cells <- list.files(file.path(cdir, "cells"), full.names = TRUE)
  # Keep only the first 3 cells (b=0, b=1, b=2).
  keep <- head(sort(cells), 3L)
  file.copy(keep, file.path(cdir2, "cells", basename(keep)))
  m_full <- utils::read.csv(file.path(cdir2, "manifest.csv"),
                            stringsAsFactors = FALSE)
  m_part <- m_full[m_full$b %in% 0:2, , drop = FALSE]
  utils::write.csv(m_part, file.path(cdir2, "manifest.csv"),
                    row.names = FALSE, quote = FALSE)

  resumed <- suppressMessages(didgpu_resume(cdir2, p))
  expect_equal(as.numeric(resumed$results$Effects[, "Estimate"]),
               as.numeric(full$results$Effects[, "Estimate"]))
  expect_equal(as.numeric(resumed$results$Effects[, "SE"]),
               as.numeric(full$results$Effects[, "SE"]))
})

test_that("didgpu_bootstrap_more extends an existing checkpoint", {
  # Run 5 reps to create a checkpoint, then add 5 more.
  # Result must match a fresh run with 10 reps and the same seed.
  p <- small_panel()
  cdir <- tempfile("didgpu_more_")
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)
  first <- didgpu(p, "Y", "unit", "period", "D",
                   effects = 2L, bootstrap_reps = 5L, seed = 11L,
                   checkpoint_dir = cdir, backend = "r", verbose = FALSE)
  extended <- suppressMessages(didgpu_bootstrap_more(cdir, p, extra_reps = 5L))

  cdir_full <- tempfile("didgpu_full_")
  on.exit(unlink(cdir_full, recursive = TRUE), add = TRUE)
  full <- didgpu(p, "Y", "unit", "period", "D",
                  effects = 2L, bootstrap_reps = 10L, seed = 11L,
                  checkpoint_dir = cdir_full, backend = "r", verbose = FALSE)

  # The extended SEs and CIs come from the same 10-rep distribution.
  expect_equal(as.numeric(extended$results$Effects[, "Estimate"]),
               as.numeric(full$results$Effects[, "Estimate"]))
  expect_equal(as.numeric(extended$results$Effects[, "SE"]),
               as.numeric(full$results$Effects[, "SE"]))
  # The meta should reflect the new larger count.
  meta <- jsonlite::fromJSON(file.path(cdir, "meta.json"),
                              simplifyVector = TRUE)
  expect_equal(meta$bootstrap_reps, 10L)
})

test_that("didgpu_resume detects panel hash mismatch", {
  p <- small_panel()
  cdir <- tempfile("didgpu_hashcheck_")
  on.exit(unlink(cdir, recursive = TRUE), add = TRUE)
  didgpu(p, "Y", "unit", "period", "D",
          effects = 1L, bootstrap_reps = 1L, seed = 1L,
          checkpoint_dir = cdir, backend = "r", verbose = FALSE)
  p2 <- p; p2$Y[1L] <- p2$Y[1L] + 1
  expect_error(didgpu_resume(cdir, p2),
               "Panel hash mismatch")
})
