// ============================================================================
// win_msvc_stubs.cpp — MinGW <-> MSVC ABI bridge for the Windows GPU build.
//
// On Windows, the .cu kernels are compiled by nvcc using MSVC (cl.exe)
// as the host compiler, but the final didgpu.dll is linked by Rtools'
// MinGW ld. MSVC emits references to a handful of runtime-support
// symbols that live in MSVC's CRT, which MinGW's runtime does not
// provide and MinGW's ld cannot pull from MSVC import libs:
//
//   __security_cookie        buffer-overrun canary (/GS)
//   __security_check_cookie  validates the canary on function exit
//   __GSHandlerCheck         SEH handler that also checks the canary
//
// We intend to compile the kernels with /GS- (no stack protector), but
// that flag does not reliably survive the MSYS make -> cmd -> nvcc
// argument mangling, so the symbols leak into the objects. Rather than
// fight the escaping, we satisfy the linker by defining no-op stubs
// here (compiled by MinGW g++ into the DLL). Functionally this is
// identical to /GS-: the buffer-overrun check is removed from the
// internal GPU helper functions, which are no-throw numeric kernels
// where it adds nothing.
//
// Guarded so the translation unit is empty everywhere except a Windows
// build with CUDA enabled (Linux, mac, and CPU-only Windows all see an
// empty file).
// ============================================================================

#if defined(_WIN32) && defined(HAS_CUDA)

extern "C" {

// Normally a per-process random value; a fixed constant is fine here
// because we never actually validate it (the check below is a no-op).
void* __security_cookie = reinterpret_cast<void*>(0x00002B992DDFA232ULL);

// MSVC calls this on function exit to compare the on-stack canary with
// __security_cookie. We skip the check entirely.
void __security_check_cookie(void* /*stack_cookie*/) {}

// Referenced from the MSVC objects' .xdata (SEH unwind) for functions
// that would carry a GS check. Never invoked in our no-throw kernels;
// the definition just satisfies the linker.
void __GSHandlerCheck(void) {}

}  // extern "C"

#endif  // _WIN32 && HAS_CUDA
