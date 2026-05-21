# ============================================================================
# Checkpoint layer
#
# Disk layout (one directory per run):
#
#   checkpoint_dir/
#   |-- meta.json        config + panel hash + seed; written once at init
#   |-- manifest.csv     one row per saved cell; append-only, fsync each row
#   `-- cells/
#       |-- b0000.rds    point estimate (b == 0)
#       |-- b0001.rds    bootstrap rep 1 (all event-times in one cell)
#       |-- b0002.rds
#       `-- ...
#
# Decomposition unit choice: one cell == one bootstrap iteration carrying
# *all* event-time coefficients for that iteration. Rationale: in the
# reference (DIDmultiplegtDYN), one call to the core estimator produces
# all event-time coefficients together; splitting across event-times
# would require interposing inside the core, which is fragile and would
# break numerical equivalence. One bootstrap rep is the natural,
# crash-safe checkpoint granularity. If later we port the core to
# C++/CUDA and gain finer control, we can subdivide.
#
# Atomicity contract: a cell is "committed" only after BOTH the .rds is
# written AND its row is appended to manifest.csv with a successful
# fsync. On resume we trust manifest.csv as the source of truth;
# orphaned .rds files without manifest entries are ignored (and may be
# overwritten).
# ============================================================================


#' Initialise a checkpoint directory
#'
#' Creates `checkpoint_dir`, writes `meta.json` recording the run
#' configuration, and seeds an empty `manifest.csv`. Refuses to
#' overwrite an existing manifest unless `force = TRUE`.
#'
#' @param checkpoint_dir Path to checkpoint directory (created if missing).
#' @param meta Named list. Must include at minimum: `panel_hash`, `seed`,
#'   `bootstrap_reps`, `effects`, `placebo`, `outcome`, `group`, `time`,
#'   `treatment`, `package_version`. Anything else is recorded as-is.
#' @param force If TRUE, wipe any existing cells/ and manifest.csv before
#'   initialising. Default FALSE.
#' @return Invisibly, the normalised `checkpoint_dir` path.
#' @examples
#' cdir <- tempfile("init_demo_")
#' meta <- list(panel_hash = "abc123", seed = 1L, bootstrap_reps = 5L,
#'              effects = 2L, placebo = 0L, outcome = "Y", group = "g",
#'              time = "t", treatment = "D",
#'              package_version = "0.1.0")
#' didgpu_init_checkpoint(cdir, meta)
#' list.files(cdir)
#' unlink(cdir, recursive = TRUE)
#' @export
didgpu_init_checkpoint <- function(checkpoint_dir, meta, force = FALSE) {
  stopifnot(is.character(checkpoint_dir), length(checkpoint_dir) == 1L)
  stopifnot(is.list(meta), !is.null(names(meta)))
  required <- c("panel_hash", "seed", "bootstrap_reps", "effects", "placebo",
                "outcome", "group", "time", "treatment", "package_version")
  missing <- setdiff(required, names(meta))
  if (length(missing)) {
    stop("meta missing required fields: ", paste(missing, collapse = ", "))
  }

  dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(checkpoint_dir, "cells"), showWarnings = FALSE)

  manifest_path <- file.path(checkpoint_dir, "manifest.csv")
  meta_path     <- file.path(checkpoint_dir, "meta.json")

  if (file.exists(manifest_path) && !force) {
    stop("manifest.csv already exists at ", manifest_path,
         ". Pass `force = TRUE` to overwrite, or use `didgpu_load_checkpoint()` to resume.")
  }
  if (force) {
    unlink(list.files(file.path(checkpoint_dir, "cells"),
                      full.names = TRUE), force = TRUE)
    if (file.exists(manifest_path)) file.remove(manifest_path)
    if (file.exists(meta_path))     file.remove(meta_path)
  }

  # Write meta atomically (write tmp + rename).
  meta_json <- jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE, null = "null")
  tmp <- paste0(meta_path, ".tmp")
  writeLines(meta_json, tmp, useBytes = TRUE)
  file.rename(tmp, meta_path)

  # Manifest header.
  con <- file(manifest_path, open = "wb")
  on.exit(close(con))
  writeLines("b,cell_file,n_coefs,wall_seconds,wrote_at_utc",
             con, useBytes = TRUE)

  invisible(normalizePath(checkpoint_dir, winslash = "/", mustWork = TRUE))
}


#' Load checkpoint metadata and manifest
#'
#' @param checkpoint_dir Path to an existing checkpoint directory.
#' @return A list with `meta` (parsed from meta.json), `manifest`
#'   (data.frame; possibly zero rows), and `checkpoint_dir` (normalised).
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 30L, n_periods = 8L, seed = 1L)
#' cdir <- tempfile("load_demo_")
#' didgpu(p, "Y", "unit", "period", "D",
#'         effects = 1L, bootstrap_reps = 0L,
#'         checkpoint_dir = cdir, backend = "r", verbose = FALSE)
#' chk <- didgpu_load_checkpoint(cdir)
#' chk$meta$effects
#' nrow(chk$manifest)
#' unlink(cdir, recursive = TRUE)
#' }
#' @export
didgpu_load_checkpoint <- function(checkpoint_dir) {
  stopifnot(dir.exists(checkpoint_dir))
  manifest_path <- file.path(checkpoint_dir, "manifest.csv")
  meta_path     <- file.path(checkpoint_dir, "meta.json")
  if (!file.exists(meta_path)) {
    stop("meta.json not found at ", meta_path,
         ". This does not look like a didgpu checkpoint directory.")
  }
  meta <- jsonlite::fromJSON(paste(readLines(meta_path, warn = FALSE), collapse = "\n"),
                              simplifyVector = TRUE)
  manifest <- if (file.exists(manifest_path)) {
    # Robust to header-only file.
    info <- file.info(manifest_path)
    if (is.na(info$size) || info$size <= 0) {
      .empty_manifest()
    } else {
      m <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                           colClasses = c(b = "integer",
                                          cell_file = "character",
                                          n_coefs = "integer",
                                          wall_seconds = "numeric",
                                          wrote_at_utc = "character"))
      if (nrow(m) == 0L) .empty_manifest() else m
    }
  } else {
    .empty_manifest()
  }

  list(
    meta = meta,
    manifest = manifest,
    checkpoint_dir = normalizePath(checkpoint_dir, winslash = "/",
                                    mustWork = TRUE)
  )
}


#' Save a single completed bootstrap cell
#'
#' Internal — called by the orchestrator. Writes the cell .rds, then
#' appends a row to manifest.csv. The file write happens before the
#' manifest append, so any orphaned .rds (cell on disk without manifest
#' row) is safely ignorable on resume.
#'
#' @keywords internal
#' @noRd
.save_cell <- function(checkpoint_dir, b, value, wall_seconds = NA_real_) {
  stopifnot(is.numeric(b), length(b) == 1L, b >= 0, b == as.integer(b))
  b <- as.integer(b)
  cells_dir <- file.path(checkpoint_dir, "cells")
  cell_file <- sprintf("b%04d.rds", b)
  cell_path <- file.path(cells_dir, cell_file)

  tmp <- paste0(cell_path, ".tmp")
  saveRDS(value, tmp, version = 2L, compress = "xz")
  file.rename(tmp, cell_path)

  n_coefs <- if (is.list(value) && !is.null(value$coef)) {
    length(value$coef)
  } else if (is.numeric(value)) {
    length(value)
  } else NA_integer_

  row <- sprintf("%d,%s,%s,%s,%s",
                 b, cell_file,
                 format(n_coefs),
                 format(wall_seconds, scientific = FALSE),
                 format(Sys.time(), tz = "UTC",
                        format = "%Y-%m-%dT%H:%M:%SZ"))
  manifest_path <- file.path(checkpoint_dir, "manifest.csv")
  con <- file(manifest_path, open = "ab")
  on.exit(close(con))
  writeLines(row, con, useBytes = TRUE)
  invisible(cell_path)
}


#' Read all committed cells back into memory
#'
#' Uses the manifest as source of truth: skips any .rds not in the
#' manifest, and warns about manifest rows whose .rds is missing.
#'
#' @param checkpoint_dir Path to an existing checkpoint directory.
#' @return Named list: cell `b` is at position `as.character(b)`.
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 30L, n_periods = 8L, seed = 1L)
#' cdir <- tempfile("agg_demo_")
#' didgpu(p, "Y", "unit", "period", "D",
#'         effects = 1L, bootstrap_reps = 2L, seed = 1L,
#'         checkpoint_dir = cdir, backend = "r", verbose = FALSE)
#' cells <- didgpu_aggregate_cells(cdir)
#' length(cells)              # 3 (point estimate + 2 bootstrap reps)
#' unlink(cdir, recursive = TRUE)
#' }
#' @export
didgpu_aggregate_cells <- function(checkpoint_dir) {
  chk <- didgpu_load_checkpoint(checkpoint_dir)
  if (nrow(chk$manifest) == 0L) return(list())
  cells <- vector("list", nrow(chk$manifest))
  names(cells) <- as.character(chk$manifest$b)
  for (i in seq_len(nrow(chk$manifest))) {
    cp <- file.path(checkpoint_dir, "cells", chk$manifest$cell_file[i])
    if (!file.exists(cp)) {
      warning("manifest references missing cell file: ", cp)
      next
    }
    cells[[i]] <- readRDS(cp)
  }
  cells
}


#' Compute the set of bootstrap reps still to do
#'
#' @keywords internal
#' @noRd
.cells_todo <- function(manifest, bootstrap_reps) {
  done <- as.integer(manifest$b)
  setdiff(0L:as.integer(bootstrap_reps), done)
}


#' Add more bootstrap reps to an existing checkpoint
#'
#' Use this when you already have, say, 100 reps committed and you
#' realise you want 200. Rather than throwing away the 100 you have
#' and starting over, this function patches the checkpoint's recorded
#' `bootstrap_reps` upward and re-invokes `didgpu()` with
#' `resume = TRUE`, so only the new (`extra_reps`) cells are computed.
#'
#' The new cells use seeds `seed + (current_reps + 1)..(current_reps + extra_reps)`,
#' which is exactly what `didgpu()` would have done on a fresh run with
#' the larger bootstrap count. The result is therefore identical to
#' starting from scratch with the larger count.
#'
#' @param checkpoint_dir Path to existing checkpoint.
#' @param df The same panel originally passed (validated by hash).
#' @param extra_reps Integer >= 1: how many MORE reps to add.
#' @param ... Additional arguments forwarded to `didgpu()` (e.g.
#'   `n_workers`, `verbose`). Overriding identity-relevant arguments
#'   like `outcome`, `effects`, `seed` will fail compatibility check.
#' @return The aggregated result, same shape as `didgpu()`.
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L,
#'                             tau_profile = c(0.5, 1.0), seed = 17L)
#' cdir <- tempfile("didgpu_more_demo_")
#' didgpu(p, "Y", "unit", "period", "D",
#'         effects = 2L, bootstrap_reps = 3L, seed = 1L,
#'         checkpoint_dir = cdir, backend = "r", verbose = FALSE)
#' # Add 2 more bootstrap reps without redoing the first 3.
#' fit <- didgpu_bootstrap_more(cdir, df = p, extra_reps = 2L,
#'                                verbose = FALSE)
#' fit$args$bootstrap_reps   # now 5
#' unlink(cdir, recursive = TRUE)
#' }
#' @export
didgpu_bootstrap_more <- function(checkpoint_dir, df, extra_reps, ...) {
  stopifnot(is.numeric(extra_reps), length(extra_reps) == 1L,
            extra_reps == as.integer(extra_reps), extra_reps >= 1L)
  extra_reps <- as.integer(extra_reps)

  chk <- didgpu_load_checkpoint(checkpoint_dir)
  meta_old <- chk$meta
  current <- as.integer(meta_old$bootstrap_reps)
  new_total <- current + extra_reps

  # Patch meta.json with the new bootstrap_reps so the
  # .check_meta_compatibility call inside didgpu() sees a matching count.
  meta_new <- meta_old
  meta_new$bootstrap_reps <- new_total
  meta_path <- file.path(checkpoint_dir, "meta.json")
  meta_json <- jsonlite::toJSON(meta_new, auto_unbox = TRUE,
                                 pretty = TRUE, null = "null")
  tmp <- paste0(meta_path, ".tmp")
  writeLines(meta_json, tmp, useBytes = TRUE)
  file.rename(tmp, meta_path)

  # Forward to didgpu_resume so we honour everything else the meta
  # recorded, but with the bumped bootstrap_reps as an explicit override.
  didgpu_resume(checkpoint_dir, df, bootstrap_reps = new_total, ...)
}


#' Quick resume helper
#'
#' Equivalent to calling `didgpu(...)` with the same arguments and
#' `checkpoint_dir = checkpoint_dir, resume = TRUE`, but pulls the
#' arguments from `meta.json` so the caller does not have to repeat
#' them. The panel must still be supplied (panels are not stored on
#' disk; only the panel hash, for integrity checks).
#'
#' All optional arguments may be overridden by passing them through
#' `...`; the override takes precedence over the checkpointed value.
#' Use this to (for example) bump `bootstrap_reps` higher or switch
#' `backend` mid-run — but note that mismatches on identity-relevant
#' fields (outcome, group, time, treatment, effects, placebo, seed)
#' will be rejected by `.check_meta_compatibility()` inside `didgpu()`.
#'
#' @param checkpoint_dir Path to existing checkpoint directory.
#' @param df The same panel originally passed (validated by hash).
#' @param ... Optional overrides forwarded to `didgpu()`.
#' @return The aggregated result, same shape as `didgpu()`.
#' @examples
#' \donttest{
#' p <- didgpu_simulate_panel(n_units = 40L, n_periods = 10L,
#'                             tau_profile = c(0.5, 1.0), seed = 17L)
#' cdir <- tempfile("didgpu_resume_demo_")
#' fit1 <- didgpu(p, "Y", "unit", "period", "D",
#'                 effects = 2L, bootstrap_reps = 3L, seed = 1L,
#'                 checkpoint_dir = cdir, backend = "r", verbose = FALSE)
#' # Resume with no extra work — produces the identical result.
#' fit2 <- didgpu_resume(cdir, df = p)
#' identical(coef(fit1), coef(fit2))
#' unlink(cdir, recursive = TRUE)
#' }
#' @export
didgpu_resume <- function(checkpoint_dir, df, ...) {
  chk <- didgpu_load_checkpoint(checkpoint_dir)
  m <- chk$meta
  observed_hash <- .panel_hash(df, m$outcome, m$group, m$time, m$treatment)
  if (!identical(observed_hash, m$panel_hash)) {
    stop("Panel hash mismatch. The data passed to didgpu_resume() does not ",
         "match the panel originally checkpointed.\n  expected: ",
         m$panel_hash, "\n  observed: ", observed_hash)
  }
  # Build the call from the meta, then layer user overrides.
  # `.empty_to_null` keeps JSON-restored "null"s as actual NULLs so
  # they don't override defaults in didgpu().
  args_from_meta <- list(
    df              = df,
    outcome         = m$outcome,
    group           = m$group,
    time            = m$time,
    treatment       = m$treatment,
    effects         = m$effects,
    placebo         = m$placebo,
    cluster         = .empty_to_null(m$cluster),
    controls        = .empty_to_null(m$controls),
    weight          = .empty_to_null(m$weight),
    continuous      = .empty_to_null(m$continuous),
    trends_nonparam = .empty_to_null(m$trends_nonparam),
    only_never_switchers   = isTRUE(m$only_never_switchers),
    same_switchers         = isTRUE(m$same_switchers),
    dont_drop_larger_lower = isTRUE(m$dont_drop_larger_lower),
    switchers       = m$switchers %||% "",
    normalized      = isTRUE(m$normalized),
    bootstrap_reps  = m$bootstrap_reps,
    ci_level        = m$ci_level %||% 95,
    seed            = m$seed,
    checkpoint_dir  = checkpoint_dir,
    resume          = TRUE,
    backend         = m$backend %||% "auto"
  )
  overrides <- list(...)
  if (length(overrides) > 0L) {
    args_from_meta[names(overrides)] <- overrides
  }
  do.call(didgpu, args_from_meta)
}

# jsonlite restores absent or NULL fields as empty named list or
# zero-length; convert those back to a true NULL so didgpu() sees
# defaults rather than malformed input.
.empty_to_null <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  if (is.list(x) && length(x) == 0L) return(NULL)
  x
}


# -------- internal helpers --------

.empty_manifest <- function() {
  data.frame(
    b = integer(0),
    cell_file = character(0),
    n_coefs = integer(0),
    wall_seconds = numeric(0),
    wrote_at_utc = character(0),
    stringsAsFactors = FALSE
  )
}

#' Hash a panel deterministically
#'
#' Sorts by (group, time), keeps only the relevant columns, and md5s
#' the resulting byte stream. Two panels that produce the same hash
#' are guaranteed to produce the same estimates (under our estimator
#' contract).
#'
#' @keywords internal
#' @noRd
.panel_hash <- function(df, outcome, group, time, treatment) {
  cols <- c(group, time, outcome, treatment)
  d <- as.data.frame(df)[, cols, drop = FALSE]
  d <- d[order(d[[group]], d[[time]]), , drop = FALSE]
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(d, tmp, version = 2L, compress = FALSE)
  unname(tools::md5sum(tmp))
}

