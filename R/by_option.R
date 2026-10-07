# ============================================================================
# didgpu(by = , by_path = ): DIDmultiplegtDYN's per-level runs.
#
# As in did_multiplegt_dyn.R (2.4.0) and did_multiplegt_by_path.R, the
# estimator runs once per level of `by` -- a group-level, time-invariant
# variable -- or once per treatment path, and the result holds one
# `by_level_b` branch per run, a combined plot and, with save_sample, one
# combined sample. Ports of did_multiplegt_by_path(), combine_plot() and
# get_colors() (MIT-licensed, Copyright (c) 2024 Diego Ciccia, Felix
# Knau, Melitine Malezieux, Doulo Sow, Clement de Chaisemartin).
# ============================================================================


# Run didgpu() once per level. `a` is didgpu()'s own argument list.
#' @keywords internal
#' @noRd
.didgpu_by_levels <- function(a) {
  df <- a$df
  by <- a$by; by_path <- a$by_path
  if (!is.null(by) && !is.null(by_path)) {
    stop("You cannot specify by and by_path options together.")
  }
  if (!is.null(by)) {
    # The reference's check (did_multiplegt_dyn_by_check): the row-average
    # of each group's SD of `by` must be 0. On a non-numeric `by` that SD
    # is NA and the reference fails outright; didgpu checks such a
    # variable by its number of distinct values per group instead.
    bv <- df[[by]]
    if (!is.numeric(bv)) bv <- as.numeric(factor(bv))
    sd_g <- stats::aggregate(bv, by = list(grp = df[[a$group]]),
                             FUN = stats::sd, na.rm = TRUE)
    names(sd_g)[2] <- "sd_by"
    m <- merge(data.frame(grp_ = df[[a$group]]), sd_g, by.x = "grp_", by.y = "grp", all.x = TRUE)
    if (!isTRUE(mean(m$sd_by, na.rm = TRUE) == 0)) {
      stop(sprintf("The variable %s specified in the by option is time-varying. That variable should be time-invariant.", by))
    }
    by_levels <- levels(factor(df[[by]]))
  } else {
    pi <- .by_path_index(a)
    df <- pi$df
    by_levels <- pi$path
    a$same_switchers <- TRUE
  }

  # Advice the reference gives once, before the levels.
  if (a$bootstrap_reps > 0L && is.null(a$continuous)) {
    message("did_multiplegt_dyn computes by default analytical standard ",
            "errors - in most cases, there is no need to use the bootstrap ",
            "option.\nBootstrapping is a much slower alternative and we ",
            "recommend it only in combination with the continuous option.")
  }
  if (a$bootstrap_reps == 0L && !is.null(a$continuous)) {
    message("You specified the continuous option without the bootstrap ",
            "option. \nPlease be aware that we recommend to compute ",
            "bootstraped standard errors when you are using the continuous ",
            "option as the analytical standard errors can be liberal in that ",
            "case.")
  }

  out <- list(args = a[setdiff(names(a), "df")], by_levels = by_levels)
  old <- options(didgpu.in_by = TRUE); on.exit(options(old), add = TRUE)
  xl_design <- list(); xl_dfs <- list()
  for (b in seq_along(by_levels)) {
    lev <- by_levels[b]
    if (!is.null(by)) {
      sub <- df[df[[by]] == lev, , drop = FALSE]
      message(sprintf("Running did_multiplegt_dyn for %s = %s", by, lev))
    } else {
      sub <- df[df$path_XX == lev |
                  (df$yet_to_switch_XX == 1 & df$baseline_XX == substr(lev, 1, 1)), ,
                drop = FALSE]
      sub$path_XX <- sub$yet_to_switch_XX <- sub$baseline_XX <- NULL
      message(sprintf("Running did_multiplegt_dyn for treatment path (%s)", lev))
    }
    la <- a
    la$df <- sub; la$by <- NULL; la$by_path <- NULL; la$graph_off <- TRUE
    if (!is.null(a$checkpoint_dir)) {
      la$checkpoint_dir <- file.path(a$checkpoint_dir,
                                     gsub("[^A-Za-z0-9_.-]", "_", paste0("level_", lev)))
    }
    fit <- do.call(didgpu, la)
    branch <- list(results = fit$results, coef = fit$coef)
    if (!is.null(fit$design)) {
      branch$design <- fit$design
      if (fit$design$design_path != "console") {
        xl_design[[paste0("Design", .by_suffix(by, lev))]] <-
          as.data.frame(unclass(fit$design$design_mat))
      }
    }
    if (!is.null(fit$date_first_switch)) {
      branch$date_first_switch <- fit$date_first_switch
      if (fit$date_first_switch$dfs_path != "console") {
        xl_dfs <- c(xl_dfs, .by_dfs_sheets(fit$date_first_switch, by, lev))
      }
    }
    if (!is.null(fit$normalized_weights)) branch$normalized_weights <- fit$normalized_weights
    branch$plot <- fit$plot
    if (!is.null(fit$save_sample)) branch$save_sample <- fit$save_sample
    branch$level <- lev
    out[[paste0("by_level_", b)]] <- branch
  }
  if (length(xl_design)) openxlsx::write.xlsx(xl_design, file = a$design[2], gridExpand = TRUE)
  if (length(xl_dfs)) openxlsx::write.xlsx(xl_dfs, file = a$date_first_switch[2], gridExpand = TRUE)

  if (isTRUE(a$save_sample)) {
    ss <- do.call(rbind, lapply(seq_along(by_levels), function(b) out[[paste0("by_level_", b)]]$save_sample))
    for (b in seq_along(by_levels)) out[[paste0("by_level_", b)]]$save_sample <- NULL
    ss <- ss[order(ss[[a$group]], ss[[a$time]]), , drop = FALSE]
    rownames(ss) <- NULL
    out$save_sample <- ss
  }
  out$plot <- .combine_plot(out)
  class(out) <- c("didgpu_result", "list")
  if (isFALSE(a$graph_off) && !is.null(out$plot)) print(out$plot)
  out
}

.by_suffix <- function(by, lev) {
  if (is.null(by)) "" else paste0(", ", abbreviate(by, 5), "=", lev)
}

.by_dfs_sheets <- function(dfs, by, lev) {
  if (dfs$dfs_opt == "by_baseline_treat") {
    stats::setNames(lapply(seq_len(dfs$levels_baseline_treat),
                           function(i) as.data.frame(unclass(dfs[[paste0("dfs_mat", i)]]))),
                    paste0("Base treat. ",
                           vapply(seq_len(dfs$levels_baseline_treat),
                                  function(i) as.character(dfs[[paste0("level", i)]]), character(1)),
                           .by_suffix(by, lev)))
  } else {
    stats::setNames(list(as.data.frame(unclass(dfs$dfs_mat))),
                    paste0("Switch. Dates", .by_suffix(by, lev)))
  }
}


# The treatment paths (did_multiplegt_by_path.R): the `by_path` most
# common paths over the first effects, from the reference's design
# table, and each row's path, baseline and not-yet-switched flag.
.by_path_index <- function(a) {
  bp <- a$by_path
  if (is.character(bp) && tolower(bp) == "all") bp <- -1
  prepped <- .prep_panel(a$df, a$outcome, a$group, a$time, a$treatment,
                         controls = a$controls, weight = a$weight,
                         trends_nonparam = a$trends_nonparam,
                         dont_drop_larger_lower = isTRUE(a$dont_drop_larger_lower),
                         continuous = a$continuous,
                         trends_lin = isTRUE(a$trends_lin),
                         drop_if_d_miss_before_first_switch =
                           isTRUE(a$drop_if_d_miss_before_first_switch))
  h <- .clamp_horizons(prepped, a$effects, a$placebo, switchers = a$switchers %||% "",
                       trends_lin = isTRUE(a$trends_lin))
  if (isTRUE(h$none_estimable)) stop(.dcdh_no_effect_msg, call. = FALSE)
  l_XX <- h$l_eff
  pan <- .dcdh_desc_panel(prepped)
  des <- .dcdh_design(pan, list(1, "console"), a$weight, l_XX)
  n_paths <- nrow(des$design_mat)
  if (bp == -1) bp <- n_paths
  if (bp > n_paths) {
    message(sprintf("You requested %.0f treatment paths, but there are only %.0f paths in your data. The program will continue with the latter number of treatment paths.",
                    bp, n_paths))
  }
  ds <- matrix(des$design_mat[1:min(bp, n_paths), ], ncol = ncol(des$design_mat),
               nrow = min(bp, n_paths))
  path <- ds[, 3]
  for (j in 1:l_XX) path <- paste0(path, ",", ds[, 3 + j])

  keep <- !is.na(pan$group) & !is.na(pan$time)
  if ("no_wt_XX" %in% names(prepped)) keep <- keep & !prepped$no_wt_XX
  pidx <- pan[keep, c("group", "time", "time_XX", "treatment_XX", "F_g_XX")]
  for (i in 0:l_XX) {
    v <- ifelse(pidx$time_XX == pidx$F_g_XX - 1 + i, pidx$treatment_XX, NA)
    agg <- stats::aggregate(v, by = list(group = pidx$group), FUN = mean, na.rm = TRUE)
    pidx[[paste0("D_fg", i)]] <- agg[[2]][match(pidx$group, agg$group)]
    pidx[[paste0("D_fg", i)]][is.nan(pidx[[paste0("D_fg", i)]])] <- NA
  }
  pidx$path_XX <- as.character(pidx$D_fg0)
  for (j in 1:l_XX) {
    dj <- pidx[[paste0("D_fg", j)]]
    pidx$path_XX <- ifelse(!is.na(dj), paste0(pidx$path_XX, ",", dj), pidx$path_XX)
  }
  pidx$yet_to_switch_XX <- as.numeric(pidx$time_XX < pidx$F_g_XX)
  pidx$baseline_XX <- substr(pidx$path_XX, 1, 1)
  pidx <- pidx[, c("group", "time", "path_XX", "yet_to_switch_XX", "baseline_XX")]
  names(pidx)[1:2] <- c(a$group, a$time)
  m <- merge(as.data.frame(a$df), pidx, by = c(a$group, a$time))
  m <- m[order(m[[a$group]], m[[a$time]]), ]
  list(df = m, path = path)
}


# Plot colours, as the reference picks them: red, blue, green, cyan,
# magenta, violet, black, orange, then random others.
.get_colors <- function(N) {
  must_color <- c(552, 26, 81, 68, 450, 640, 24, 498)
  other_color <- setdiff(1:657, must_color)
  idx <- if (N > length(must_color)) c(must_color, sample(other_color, N - length(must_color)))
         else must_color[1:N]
  grDevices::colors()[idx]
}

# The combined plot of a by / by_path result (combine_plot()).
.combine_plot <- function(obj) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("cowplot", quietly = TRUE)) return(NULL)
  levs <- obj$by_levels
  if (!is.null(obj$args[["by"]])) {
    cols <- .get_colors(length(levs))
    lab <- paste0(obj$args$by, " = ", levs)
    p <- ggplot2::ggplot()
    for (j in seq_along(levs)) {
      dj <- obj[[paste0("by_level_", j)]]$plot$data
      if (is.null(dj)) next
      dj$lab_XX <- lab[j]
      p <- p +
        ggplot2::geom_point(data = dj, ggplot2::aes(x = Time, y = Estimate), colour = cols[j]) +
        ggplot2::geom_errorbar(data = dj, ggplot2::aes(x = Time, ymin = LB.CI, ymax = UB.CI),
                               position = ggplot2::position_dodge(0.05), width = 0.2,
                               colour = cols[j]) +
        ggplot2::geom_line(data = dj, ggplot2::aes(x = Time, y = Estimate, colour = lab_XX))
    }
    p + ggplot2::ylab("Estimate") +
      ggplot2::ggtitle("DID, from last period before treatment changes (t=0) to t") +
      ggplot2::xlab("Relative time to last period before treatment changes (t=0)") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "bottom") +
      ggplot2::scale_colour_manual("", breaks = lab, values = cols)
  } else {
    if (length(levs) > 100) {
      message("The command allows a maximum of 100 graphs to be combined in a 10 x 10 window. The resulting graph will be restricted to the first 100 treatment paths.")
    }
    sides <- ceiling(sqrt(length(levs))); len <- 1 / sides
    p <- cowplot::ggdraw()
    for (j in seq_along(levs)) {
      br <- obj[[paste0("by_level_", j)]]
      if (is.null(br$plot)) next
      p <- p + cowplot::draw_plot(
        br$plot + ggplot2::ggtitle(sprintf("Treatment path (%s); %.0f switchers", levs[j],
                                           br$results$Effects[1, 6])) + ggplot2::xlab(" "),
        width = len, height = len, y = (sides - ceiling(j / sides)) * len,
        x = ((j - 1) %% sides) * len)
    }
    p + ggplot2::ggtitle("DID from last period before treatment changes (t = 0) to t") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
}
