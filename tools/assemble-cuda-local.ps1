# Assemble a user-local CUDA 12.6 toolkit from NVIDIA redist archives.
# No admin: download component ZIPs + merge into C:\Users\ammonsj\cuda_local.
# Uses robocopy for the merge (retries on transient AV/Defender locks).
$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'
$base = "https://developer.download.nvidia.com/compute/cuda/redist/"
$dl   = "C:\Users\ammonsj\cuda_dl"
$root = "C:\Users\ammonsj\cuda_local"

# Clean prior partial merge + extraction temps; KEEP downloaded zips.
if (Test-Path $root) { Remove-Item -Recurse -Force $root }
New-Item -ItemType Directory -Force -Path $dl, $root | Out-Null
Get-ChildItem $dl -Directory -Filter "x_*" -EA SilentlyContinue | Remove-Item -Recurse -Force

$paths = @(
  "cuda_nvcc/windows-x86_64/cuda_nvcc-windows-x86_64-12.6.85-archive.zip",
  "cuda_cudart/windows-x86_64/cuda_cudart-windows-x86_64-12.6.77-archive.zip",
  "cuda_cccl/windows-x86_64/cuda_cccl-windows-x86_64-12.6.77-archive.zip",
  "libcublas/windows-x86_64/libcublas-windows-x86_64-12.6.4.1-archive.zip",
  "libcusolver/windows-x86_64/libcusolver-windows-x86_64-11.7.1.2-archive.zip",
  "libcurand/windows-x86_64/libcurand-windows-x86_64-10.3.7.77-archive.zip",
  # cuSPARSE is a RUNTIME dependency of cuSOLVER (cusolver64_11.dll imports
  # cusparse64_12.dll), and cuSPARSE in turn pulls nvJitLink. Without these
  # two, didgpu_cuda.dll links fine but fails to LOAD at runtime with
  # "LoadLibrary failure: The specified module could not be found".
  "libcusparse/windows-x86_64/libcusparse-windows-x86_64-12.5.4.2-archive.zip",
  "libnvjitlink/windows-x86_64/libnvjitlink-windows-x86_64-12.6.85-archive.zip"
)

foreach ($p in $paths) {
  $fn   = Split-Path $p -Leaf
  $dest = Join-Path $dl $fn
  if ((-not (Test-Path $dest)) -or ((Get-Item $dest).Length -lt 1024)) {
    Write-Output "Downloading $fn ..."
    Invoke-WebRequest -Uri "$base$p" -OutFile $dest -TimeoutSec 1800
  } else { Write-Output "Have $fn (skip download)" }
  $tmp = Join-Path $dl ("x_" + [IO.Path]::GetFileNameWithoutExtension($fn))
  if (Test-Path $tmp) { Remove-Item -Recurse -Force $tmp }
  Write-Output "Extracting $fn ..."
  Expand-Archive -Path $dest -DestinationPath $tmp -Force
  $inner = (Get-ChildItem $tmp -Directory | Select-Object -First 1).FullName
  Write-Output "Merging into cuda_local (robocopy) ..."
  # robocopy: /E recurse incl empty, quiet flags, retry 5x wait 2s.
  $null = robocopy $inner $root /E /NFL /NDL /NJH /NJS /R:5 /W:2
  if ($LASTEXITCODE -ge 8) { throw "robocopy failed for $fn (code $LASTEXITCODE)" }
}

Write-Output "`n=== Verify assembled tree ==="
"nvcc.exe:        " + (Test-Path "$root\bin\nvcc.exe")
"crt/host_config: " + (Test-Path "$root\include\crt\host_config.h")
"curand.lib:      " + (Test-Path "$root\lib\x64\curand.lib")
"cublas.lib:      " + (Test-Path "$root\lib\x64\cublas.lib")
"cusolver.lib:    " + (Test-Path "$root\lib\x64\cusolver.lib")
"cudart.lib:      " + (Test-Path "$root\lib\x64\cudart.lib")
"curand.h:        " + (Test-Path "$root\include\curand.h")
"cudart64 DLL:    " + ((Get-ChildItem "$root\bin\cudart64_*.dll" -EA SilentlyContinue | Measure-Object).Count)
"curand64 DLL:    " + ((Get-ChildItem "$root\bin\curand64_*.dll" -EA SilentlyContinue | Measure-Object).Count)
"cublas64 DLL:    " + ((Get-ChildItem "$root\bin\cublas64_*.dll" -EA SilentlyContinue | Measure-Object).Count)
"cusolver64 DLL:  " + ((Get-ChildItem "$root\bin\cusolver64_*.dll" -EA SilentlyContinue | Measure-Object).Count)
"cusparse64 DLL:  " + ((Get-ChildItem "$root\bin\cusparse64_*.dll" -EA SilentlyContinue | Measure-Object).Count)
"nvJitLink DLL:   " + ((Get-ChildItem "$root\bin\nvJitLink*.dll" -EA SilentlyContinue | Measure-Object).Count)
Write-Output "DONE"
