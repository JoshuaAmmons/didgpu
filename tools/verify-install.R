library(didgpu)
cat("loaded: OK\n")
cat("CUDA compiled in: ", didgpu_has_cuda_support(), "\n", sep = "")
if (didgpu_has_cuda_support()) {
  r <- didgpu_run_saxpy(2.0, as.numeric(1:5), as.numeric(1:5))
  cat("SAXPY result: ", paste(r, collapse = ", "), "\n", sep = "")
  cat("SAXPY expected: 3, 6, 9, 12, 15\n")
  stopifnot(identical(r, c(3, 6, 9, 12, 15)))
  cat("SAXPY: PASS\n")
}
