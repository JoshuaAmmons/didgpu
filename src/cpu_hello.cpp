// Smoke target for the CPU build path. Returns a fixed value so we can
// verify that src/cpu compilation works in isolation from any GPU
// machinery. Replace once the real per-cell solve is ported.
#include <Rcpp.h>

// [[Rcpp::export]]
double didgpu_cpu_hello() {
  return 42.0;
}
