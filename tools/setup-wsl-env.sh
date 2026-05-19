#!/usr/bin/env bash
# ============================================================================
# setup-wsl-env.sh
#
# One-shot installer for the WSL2 Ubuntu dev environment for didgpu.
# Idempotent — safe to re-run; skips anything already installed.
#
# What it installs:
#   1) Build tools (build-essential, make, git)
#   2) CUDA Toolkit 12.6 (WSL-Ubuntu variant — uses Windows-side driver)
#   3) R 4.4.x from CRAN's Ubuntu repo
#   4) R-dev support libraries needed by Rcpp/RcppEigen
#   5) R packages: Rcpp, RcppEigen, data.table, testthat, broom
#
# Usage (from PowerShell or cmd):
#   wsl.exe -d Ubuntu -- bash "/mnt/c/Users/ammonsj/DID GPU/didgpu/tools/setup-wsl-env.sh"
#
# You will be prompted for your Ubuntu sudo password once.
# Total runtime: ~15-25 min depending on network speed.
# Total download: ~2-3 GB (most of that is CUDA Toolkit).
# ============================================================================

set -euo pipefail

# Color output helpers.
red()   { printf "\033[31m%s\033[0m\n" "$*"; }
green() { printf "\033[32m%s\033[0m\n" "$*"; }
blue()  { printf "\033[34m%s\033[0m\n" "$*"; }

step() { blue ""; blue "=== $* ==="; }

# ----------------------------------------------------------------------------
# 0) Pre-flight checks
# ----------------------------------------------------------------------------
step "Pre-flight checks"

if ! command -v wsl.exe >/dev/null 2>&1; then
  # Probably means we're running inside WSL already — good.
  :
fi

# Confirm we're on Ubuntu and have nvidia-smi (driver passthrough).
. /etc/os-release
green "OS: $PRETTY_NAME"

if ! command -v nvidia-smi >/dev/null 2>&1; then
  # On WSL, nvidia-smi may live at /usr/lib/wsl/lib/nvidia-smi.
  if [ -x /usr/lib/wsl/lib/nvidia-smi ]; then
    green "nvidia-smi found at /usr/lib/wsl/lib/nvidia-smi"
  else
    red "nvidia-smi not found — WSL2 GPU passthrough may be broken."
    red "Make sure the Windows NVIDIA driver is recent (>=535) and WSL is v2."
    exit 1
  fi
else
  green "nvidia-smi: $(nvidia-smi --version 2>&1 | head -1)"
fi

# ----------------------------------------------------------------------------
# 1) Base build tools
# ----------------------------------------------------------------------------
step "Installing base build tools (build-essential, make, git, curl, wget)"

sudo apt-get update -qq
sudo apt-get install -y --no-install-recommends \
  build-essential \
  make \
  git \
  curl \
  wget \
  ca-certificates \
  software-properties-common \
  dirmngr \
  gnupg \
  pkg-config

green "Base build tools installed."

# ----------------------------------------------------------------------------
# 2) CUDA Toolkit 12.6 (WSL-Ubuntu variant)
# ----------------------------------------------------------------------------
step "Installing CUDA Toolkit 12.6 (WSL-Ubuntu variant)"

# The WSL-Ubuntu variant does NOT include the driver — the driver lives on
# the Windows side and is exposed to WSL via /usr/lib/wsl/lib/. This avoids
# the driver-version conflicts that plague generic Linux CUDA installs in WSL.

if [ ! -f /usr/local/cuda/bin/nvcc ] && ! command -v nvcc >/dev/null 2>&1; then
  TMPDIR=$(mktemp -d)
  trap 'rm -rf "$TMPDIR"' EXIT

  cd "$TMPDIR"
  wget -q https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
  sudo dpkg -i cuda-keyring_1.1-1_all.deb
  sudo apt-get update -qq
  sudo apt-get install -y --no-install-recommends cuda-toolkit-12-6

  green "CUDA Toolkit 12.6 installed at /usr/local/cuda-12.6"
else
  green "nvcc already present — skipping CUDA install."
fi

# Add CUDA to PATH for this shell session and persistently via ~/.bashrc.
if ! grep -q "/usr/local/cuda/bin" ~/.bashrc 2>/dev/null; then
  cat >> ~/.bashrc <<'EOF'

# Added by didgpu/tools/setup-wsl-env.sh — CUDA Toolkit on PATH.
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}
EOF
  green "Added CUDA to PATH in ~/.bashrc."
fi
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}

# Verify nvcc works.
if ! nvcc --version >/dev/null 2>&1; then
  red "nvcc not on PATH after install — manual fix needed."
  exit 1
fi
green "nvcc: $(nvcc --version | tail -1)"

# ----------------------------------------------------------------------------
# 3) R from CRAN's Ubuntu repo
# ----------------------------------------------------------------------------
step "Installing R from CRAN's Ubuntu repo"

if ! command -v R >/dev/null 2>&1; then
  # Add CRAN GPG key.
  if [ ! -f /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc ]; then
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
      | sudo tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc >/dev/null
  fi

  # Add CRAN apt repo for the running Ubuntu release.
  UBUNTU_CODENAME=$(lsb_release -cs)
  if ! grep -q "cloud.r-project.org" /etc/apt/sources.list.d/*.list 2>/dev/null; then
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${UBUNTU_CODENAME}-cran40/" \
      | sudo tee /etc/apt/sources.list.d/cran.list >/dev/null
  fi

  sudo apt-get update -qq
  sudo apt-get install -y --no-install-recommends r-base r-base-dev

  green "R installed: $(R --version | head -1)"
else
  green "R already present: $(R --version | head -1)"
fi

# Headers required by common R packages (Rcpp, curl, openssl, xml2).
sudo apt-get install -y --no-install-recommends \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libfontconfig1-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev

green "R-dev support libraries installed."

# ----------------------------------------------------------------------------
# 4) R packages didgpu needs
# ----------------------------------------------------------------------------
step "Installing R packages (Rcpp, RcppEigen, data.table, testthat, broom)"

# Install into the system-wide site-library so any user (including root)
# can load them. This requires sudo. The alternative is a per-user library,
# which works but is fragile — you'd lose package access when running as
# a different WSL user or via wsl -u root.
sudo R --no-save <<'EOF' || { echo "R package install failed"; exit 1; }
options(repos = c(CRAN = "https://cloud.r-project.org"))
needed <- c("Rcpp", "RcppEigen", "data.table", "testthat", "broom",
            "quadprog", "ggplot2")
to_install <- setdiff(needed, rownames(installed.packages()))
if (length(to_install)) {
  install.packages(to_install, Ncpus = parallel::detectCores())
} else {
  cat("All R packages already installed.\n")
}
EOF

green "R packages installed."

# ----------------------------------------------------------------------------
# Done.
# ----------------------------------------------------------------------------
step "Setup complete"
green ""
green "Environment:"
green "  $(gcc --version | head -1)"
green "  $(g++ --version | head -1)"
green "  $(nvcc --version | tail -1)"
green "  $(R --version | head -1)"
green ""
green "GPU (from WSL):"
nvidia-smi -L 2>/dev/null || /usr/lib/wsl/lib/nvidia-smi -L
green ""
green "Next step: run from a fresh WSL shell so PATH picks up CUDA:"
green "    wsl.exe -d Ubuntu"
green "    cd '/mnt/c/Users/ammonsj/DID GPU/didgpu'"
green "    R CMD INSTALL ."
green ""
