#!/usr/bin/env bash
# ============================================================================
# finish-r-install.sh
#
# Cleanup after setup-wsl-env.sh ran but 3 R packages failed because
# libuv1-dev was missing. Installs libuv1-dev, then retries the failed
# R packages (fs, pkgload, testthat — and anything else that needs them).
#
# Idempotent — safe to re-run.
# ============================================================================

set -euo pipefail

echo "=== Installing libuv1-dev ==="
sudo apt-get update -qq
sudo apt-get install -y --no-install-recommends libuv1-dev

echo ""
echo "=== Retrying R packages ==="
sudo R --no-save <<'EOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))
needed <- c("fs", "pkgload", "testthat",
            "Rcpp", "RcppEigen", "data.table", "broom",
            "quadprog", "ggplot2")
to_install <- setdiff(needed, rownames(installed.packages()))
if (length(to_install)) {
  cat("Installing:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, Ncpus = parallel::detectCores())
} else {
  cat("All packages already installed.\n")
}
# Final sanity check:
still_missing <- setdiff(needed, rownames(installed.packages()))
if (length(still_missing)) {
  stop("Still missing: ", paste(still_missing, collapse = ", "))
}
cat("\nAll required R packages installed:\n")
print(installed.packages()[needed, "Version", drop = FALSE])
EOF

echo ""
echo "=== Done ==="
