#!/usr/bin/env bash
# ============================================================================
# build-on-wsl.sh
#
# Build the didgpu package inside WSL2 Ubuntu. Phase 0b of the roadmap:
# verify the package compiles end-to-end with the CUDA backend enabled.
#
# Run after setup-wsl-env.sh has installed gcc, nvcc, R, and the R deps.
#
# Usage from PowerShell:
#   wsl.exe -d Ubuntu -- bash "/mnt/c/Users/ammonsj/DID GPU/didgpu/tools/build-on-wsl.sh"
# ============================================================================

set -euo pipefail

# Color helpers.
green() { printf "\033[32m%s\033[0m\n" "$*"; }
blue()  { printf "\033[34m%s\033[0m\n" "$*"; }
red()   { printf "\033[31m%s\033[0m\n" "$*"; }
step()  { blue ""; blue "=== $* ==="; }

PKG_DIR="/mnt/c/Users/ammonsj/DID GPU/didgpu"
cd "$PKG_DIR"

# Make sure CUDA is on PATH for this session (setup-wsl-env.sh adds it
# to ~/.bashrc, but a non-interactive bash invocation may not source it).
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}

step "Toolchain check"
echo "  gcc:  $(gcc --version | head -1)"
echo "  nvcc: $(nvcc --version | tail -1)"
echo "  R:    $(R --version | head -1)"
echo "  GPU:  $(nvidia-smi -L 2>/dev/null || /usr/lib/wsl/lib/nvidia-smi -L)"

step "Cleaning previous build artifacts in src/"
rm -f src/*.o src/*.so src/symbols.rds 2>/dev/null || true

step "R CMD INSTALL ."
# We install as root because the package R libs live in
# /usr/local/lib/R/site-library on Ubuntu (the standard CRAN package
# location), and only root can write there. Same pattern as
# finish-r-install.sh.
sudo R CMD INSTALL --no-multiarch . 2>&1 | tee /tmp/didgpu-install.log

step "Smoke test"
sudo R --no-save <<'EOF'
library(didgpu)
cat("Loaded didgpu\n")
cat("CUDA support compiled in?", didgpu_has_cuda_support(), "\n")
# If CUDA is wired through, run a tiny SAXPY end-to-end.
if (didgpu_has_cuda_support()) {
  cat("Trying didgpu_run_saxpy(2.0, x = 1:5, y = 1:5)...\n")
  result <- tryCatch(
    didgpu_run_saxpy(2.0, x = as.numeric(1:5), y = as.numeric(1:5)),
    error = function(e) paste("ERROR:", conditionMessage(e))
  )
  cat("Result:", paste(result, collapse = ", "), "\n")
  cat("Expected: 3, 6, 9, 12, 15 (since 2*x + y for x=y=1:5)\n")
}
# Always run a tiny R-backend smoke test.
cat("\nR backend smoke test:\n")
p <- didgpu_simulate_panel(n_units = 20L, n_periods = 6L,
                           tau_profile = c(0.5, 1.0), seed = 17L)
fit <- didgpu(p, "Y", "unit", "period", "D",
              effects = 2L, bootstrap_reps = 0L,
              backend = "r", verbose = FALSE)
cat("ATE:", fit$results$ATE[1, "Estimate"], "\n")
cat("Sample size:", attr(fit, "n_rows"), "\n")
EOF

green ""
green "=== Build complete ==="
green "Install log: /tmp/didgpu-install.log"
