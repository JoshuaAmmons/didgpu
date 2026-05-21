@echo off
REM ===========================================================================
REM nvcc_wrapper.bat
REM
REM Bridges GNU Make (invoked by R CMD INSTALL via Rtools) to nvcc on Windows.
REM nvcc requires MSVC cl.exe to be on PATH AND requires MSVC's INCLUDE/LIB
REM environment to be set up. vcvars64.bat does both. This wrapper:
REM   1) Locates vcvars64.bat by probing standard VS install paths.
REM   2) Calls it (silently).
REM   3) Forwards all arguments to nvcc.exe.
REM
REM First argument MUST be the full path to nvcc.exe; the rest are nvcc args.
REM ===========================================================================

setlocal EnableDelayedExpansion

REM --- find vcvars64.bat ------------------------------------------------------
REM `for /d` with spaced wildcards is buggy in cmd, so just probe an explicit
REM list of standard install locations. Covers 2017/2019/2022 x BuildTools/
REM Community/Professional/Enterprise x x86/x64 install root.
set "VCVARS="
call :probe "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
call :probe "C:\Program Files\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
call :probe "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call :probe "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat"
call :probe "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
call :probe "C:\Program Files (x86)\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call :probe "C:\Program Files (x86)\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat"
call :probe "C:\Program Files (x86)\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat"
call :probe "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
call :probe "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
call :probe "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvars64.bat"

if "!VCVARS!"=="" (
    echo nvcc_wrapper.bat: cannot locate vcvars64.bat 1>&2
    echo nvcc_wrapper.bat: install Visual Studio Build Tools 2022 with the 1>&2
    echo nvcc_wrapper.bat: "Desktop development with C++" workload. 1>&2
    exit /b 1
)

REM --- call vcvars then forward to nvcc --------------------------------------
echo nvcc_wrapper.bat: VCVARS=!VCVARS! 1>&2
call "!VCVARS!" 1>&2
if errorlevel 1 (
    echo nvcc_wrapper.bat: vcvars64.bat returned errorlevel %errorlevel% 1>&2
    exit /b 1
)
echo nvcc_wrapper.bat: about to run nvcc with: %* 1>&2

REM Run nvcc with all forwarded args.
%*
set "RC=%errorlevel%"
echo nvcc_wrapper.bat: nvcc returned %RC% 1>&2
exit /b %RC%

REM ---------------------------------------------------------------------------
:probe
REM Sets VCVARS if not already set AND %~1 exists.
if not "!VCVARS!"=="" goto :eof
if exist %1 set "VCVARS=%~1"
goto :eof
