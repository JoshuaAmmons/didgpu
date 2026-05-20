# Build didgpu on Windows R with the user-local CUDA 12.6 toolkit.
# No admin: Rtools44 + VS Build Tools + cuda_local are all user-space.
$ErrorActionPreference = 'Continue'
$repo = "C:\Users\ammonsj\DID GPU\didgpu"
$rlib = "C:\Users\ammonsj\AppData\Local\R\win-library\4.4"
$R    = "C:\Program Files\R\R-4.4.1\bin\R.exe"

# Point the build at the user-local toolkit (Windows path; Makevars.win
# cygpath-converts it). Rtools44 on PATH for make + MinGW g++.
$env:CUDA_PATH = "C:\Users\ammonsj\cuda_local"
$env:PATH = "C:\rtools44\usr\bin;C:\rtools44\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH

Set-Location $repo
Write-Output "Cleaning stale build artifacts in src/ ..."
Remove-Item "$repo\src\*.o","$repo\src\*.so","$repo\src\*.a","$repo\src\*.dll" -Force -ErrorAction SilentlyContinue

Write-Output "CUDA_PATH = $env:CUDA_PATH"
Write-Output "Running R CMD INSTALL (this compiles the .cu kernels with nvcc) ...`n"
& $R CMD INSTALL --no-multiarch --library="$rlib" "$repo" 2>&1
Write-Output "`nR CMD INSTALL exit code: $LASTEXITCODE"
