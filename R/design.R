# ============================================================================
# design and date_first_switch: DIDmultiplegtDYN's descriptive tables of
# the switchers' treatment paths and of their switching dates.
#
# Ports of did_multiplegt_dyn_design() and did_multiplegt_dyn_dfs() from
# DIDmultiplegtDYN 2.4.0 (MIT-licensed, Copyright (c) 2024 Diego Ciccia,
# Felix Knau, Melitine Malezieux, Doulo Sow, Clement de Chaisemartin).
# Both read the estimation panel the reference holds at the end of its
# run; `panel` is didgpu's equivalent, built by .dcdh_desc_panel().
# ============================================================================


# The columns both tables read, from the prepared panel: indices, the
# user's labels, switch date, baseline, treatment and the user's weight
# (NA on the balancing fill-in rows, as in the reference).
#' @keywords internal
#' @noRd
.dcdh_desc_panel <- function(prepped) {
  d <- prepped
  data.frame(group_XX = d$group_XX, time_XX = d$time_XX,
             group = d$grp_orig_XX, time = d$time_orig_XX,
             F_g_XX = d$F_g_XX, d_sq_XX = d$d_sq_XX,
             treatment_XX = d$treatment_XX,
             weight_XX = if ("wt_in_XX" %in% names(d)) d$wt_in_XX else 1,
             stringsAsFactors = FALSE)
}


#' @keywords internal
#' @noRd
.dcdh_design <- function(df, design_opt, weight, l_XX) {
  suppressWarnings({
    des_p <- as.numeric(design_opt[1])
    des_path <- design_opt[2]
    des_n <- l_XX
    des_per <- des_p * 100
    df$F_g_plus_n_XX <- df$F_g_XX + des_n - 1
    df <- subset(df, df$time_XX >= df$F_g_XX - 1 & df$time_XX <= df$F_g_plus_n_XX)
    df <- df[order(df$group_XX, df$time_XX), ]
    df$time_l_XX <- stats::ave(seq_len(nrow(df)), df$group_XX, FUN = seq_along)
    df <- subset(df, select = c("group_XX", "time_l_XX", "weight_XX",
                                "treatment_XX", "F_g_XX"))
    if (!is.null(weight)) {
      weight_sum <- stats::aggregate(df$weight_XX, by = list(group_XX = df$group_XX),
                                     FUN = sum, na.rm = TRUE)
      names(weight_sum)[2] <- "g_weight_XX"
      df <- merge(df, weight_sum, by = "group_XX", all.x = TRUE)
    } else {
      df$g_weight_XX <- 1
    }
    df$weight_XX <- NULL
    max_time <- max(df$time_l_XX, na.rm = TRUE)
    treat_list <- c()
    for (i in 1:max_time) {
      df_sub <- df[df$time_l_XX == i, c("group_XX", "treatment_XX")]
      treat_mean <- stats::aggregate(df_sub$treatment_XX,
                                     by = list(group_XX = df_sub$group_XX),
                                     FUN = mean, na.rm = TRUE)
      names(treat_mean)[2] <- paste0("treatment_XX", i)
      df <- merge(df, treat_mean, by = "group_XX", all.x = TRUE)
      treat_list <- c(treat_list, paste0("treatment_XX", i))
    }
    df$time_l_XX <- df$treatment_XX <- NULL
    df <- unique(df)
    for (var in treat_list) df <- subset(df, !is.na(df[[var]]))
    df$N_XX <- 1
    df$N_w_XX <- (df$g_weight_XX * df$N_XX) / sum(df$g_weight_XX, na.rm = TRUE)
    df$group_XX <- df$g_weight_XX <- NULL
    N_sum <- stats::aggregate(df$N_XX, by = df[treat_list], FUN = sum, na.rm = TRUE)
    names(N_sum)[ncol(N_sum)] <- "N_XX_sum"
    Nw_sum <- stats::aggregate(df$N_w_XX, by = df[treat_list], FUN = sum, na.rm = TRUE)
    names(Nw_sum)[ncol(Nw_sum)] <- "N_w_XX_sum"
    df <- merge(df, N_sum, by = treat_list, all.x = TRUE)
    df <- merge(df, Nw_sum, by = treat_list, all.x = TRUE)
    df$N_XX <- df$N_XX_sum
    df$N_w_XX <- df$N_w_XX_sum
    df$N_XX_sum <- df$N_w_XX_sum <- NULL
    df$F_g_XX <- NULL
    df <- unique(df)
    tot_switch <- sum(df$N_XX, na.rm = TRUE)
    df$neg_N_XX <- -df$N_XX
    df$treat_GRP <- as.numeric(factor(do.call(paste, c(df[treat_list], sep = "_"))))
    df <- df[order(df$neg_N_XX, df$treat_GRP), ]
    df$neg_N_XX <- df$treat_GRP <- NULL
    df$cum_sum_XX <- cumsum(df$N_w_XX)
    df$in_table_XX <- as.numeric(df$cum_sum_XX <= des_p)
    df <- df[order(df$in_table_XX, df$cum_sum_XX), ]
    df$id_XX <- stats::ave(seq_len(nrow(df)), df$in_table_XX, FUN = seq_along)
    df <- subset(df, df$in_table_XX == 1 | (df$in_table_XX == 0 & df$id_XX == 1))
    last_p <- if (des_p < 1) 100 * min(df$cum_sum_XX[df$in_table_XX == 0]) else 100
    df$neg_N_XX <- -df$N_XX
    df$treat_GRP <- as.numeric(factor(do.call(paste, c(df[treat_list], sep = "_"))))
    df <- df[order(df$neg_N_XX, df$treat_GRP), ]
    df <- subset(df, select = c("N_XX", "N_w_XX", treat_list))
    df$N_w_XX <- df$N_w_XX * 100
    coln <- c("N", "Share")
    rown <- c()
    desmat <- matrix(NA, nrow = dim(df)[1], ncol = 2 + 1 + l_XX)
    df <- data.frame(df)
    for (j in 1:(2 + 1 + l_XX)) {
      for (i in 1:dim(df)[1]) {
        if (j == 1) rown <- c(rown, paste0("TreatPath", i))
        desmat[i, j] <- as.numeric(df[i, j])
      }
      if (j > 2) coln <- c(coln, paste0("ℓ", "=", j - 2 - 1))
    }
    colnames(desmat) <- coln
    rownames(desmat) <- rown
    desmat[, 2] <- noquote(sprintf("%s", format(round(desmat[, 2], 2), big.mark = ",",
                                                scientific = FALSE, trim = TRUE)))
    des_const <- c(l_XX, des_per, tot_switch, last_p)
    names(des_const) <- c("effects", "coverage_opt", "switchers", "detected_coverage")
    list(design_path = des_path, design_mat = noquote(desmat),
         design_const = des_const)
  })
}


#' @keywords internal
#' @noRd
.dcdh_dfs <- function(df, dfs, T_max_XX) {
  dfs_opt <- dfs[1]
  dfs_path <- dfs[2]
  if (dfs_opt != "" & dfs_opt != "by_baseline_treat") {
    stop("Only option by_baseline_treat allowed.")
  }
  suppressWarnings({
    df <- subset(df, !(df$F_g_XX == T_max_XX + 1 | is.na(df$F_g_XX)))
    df <- subset(df, df$time_XX == df$F_g_XX)
    df <- subset(df, select = c("group", "time", "F_g_XX", "d_sq_XX"))
    .tab <- function(dd) {
      dfsmat <- matrix(NA, ncol = 2, nrow = dim(dd)[1])
      rown <- c()
      dd <- data.frame(dd)
      for (j in 1:2) {
        for (i in 1:dim(dd)[1]) {
          if (j == 1) rown <- c(rown, dd$time[i])
          dfsmat[i, j] <- as.numeric(dd[i, j])
        }
      }
      colnames(dfsmat) <- c("N", "Share")
      rownames(dfsmat) <- rown
      dfsmat[, 2] <- sprintf("%s", format(round(dfsmat[, 2], 2), big.mark = ",",
                                          scientific = FALSE, trim = TRUE))
      noquote(dfsmat)
    }
    if (dfs_opt == "") {
      df$tot_s <- 1
      tot_agg <- stats::aggregate(df$tot_s, by = list(time = df$time), FUN = sum, na.rm = TRUE)
      names(tot_agg)[2] <- "tot_s_sum"
      df <- merge(df, tot_agg, by = "time", all.x = TRUE)
      df$tot_s <- df$tot_s_sum
      df$tot_s_sum <- NULL
      df$group <- df$F_g_XX <- df$d_sq_XX <- NULL
      df <- unique(df)
      df <- df[order(df$time), ]
      df$share_XX <- (df$tot_s / sum(df$tot_s, na.rm = TRUE)) * 100
      df <- subset(df, select = c("tot_s", "share_XX", "time"))
      res_dfs <- list(dfs_opt = dfs_opt, dfs_path = dfs_path, dfs_mat = .tab(df))
    } else {
      df$tot_s <- 1
      tot_agg2 <- stats::aggregate(df$tot_s, by = list(time = df$time, d_sq_XX = df$d_sq_XX),
                                   FUN = sum, na.rm = TRUE)
      names(tot_agg2)[3] <- "tot_s_sum"
      df <- merge(df, tot_agg2, by = c("time", "d_sq_XX"), all.x = TRUE)
      df$tot_s <- df$tot_s_sum
      df$tot_s_sum <- NULL
      df$group <- df$F_g_XX <- NULL
      df <- unique(df)
      df <- df[order(df$d_sq_XX, df$time), ]
      levels_d_sq_XX <- levels(factor(df$d_sq_XX))
      res_dfs <- list(dfs_opt = dfs_opt, dfs_path = dfs_path,
                      levels_baseline_treat = length(levels_d_sq_XX))
      index <- 1
      for (l in levels_d_sq_XX) {
        df_by <- subset(df, df$d_sq_XX == l)
        df_by$share_XX <- (df_by$tot_s / sum(df_by$tot_s, na.rm = TRUE)) * 100
        df_by <- subset(df_by, select = c("tot_s", "share_XX", "time"))
        res_dfs <- append(res_dfs, l)
        res_dfs <- append(res_dfs, list(.tab(df_by)))
        names(res_dfs)[length(res_dfs) - 1] <- paste0("level", index)
        names(res_dfs)[length(res_dfs)] <- paste0("dfs_mat", index)
        index <- index + 1
      }
    }
  })
  res_dfs
}


# Write the tables to the Excel file the user named, with the
# reference's sheet names (did_multiplegt_dyn.R:334-339).
#' @keywords internal
#' @noRd
.dcdh_write_xlsx <- function(design, dfs, design_opt, dfs_opt) {
  if ((!is.null(design) && design$design_path != "console") ||
      (!is.null(dfs) && dfs$dfs_path != "console")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Writing the design / date_first_switch tables to a file needs ",
           "the openxlsx package.")
    }
  }
  if (!is.null(design) && design$design_path != "console") {
    openxlsx::write.xlsx(list(Design = as.data.frame(unclass(design$design_mat))),
                         file = design_opt[2], gridExpand = TRUE)
  }
  if (!is.null(dfs) && dfs$dfs_path != "console") {
    sheets <- if (dfs$dfs_opt == "by_baseline_treat") {
      stats::setNames(lapply(seq_len(dfs$levels_baseline_treat),
                             function(i) as.data.frame(unclass(dfs[[paste0("dfs_mat", i)]]))),
                      paste0("Base treat. ", vapply(seq_len(dfs$levels_baseline_treat),
                                                    function(i) as.character(dfs[[paste0("level", i)]]),
                                                    character(1))))
    } else list(`Switch. Dates` = as.data.frame(unclass(dfs$dfs_mat)))
    openxlsx::write.xlsx(sheets, file = dfs_opt[2], gridExpand = TRUE)
  }
  invisible(NULL)
}
