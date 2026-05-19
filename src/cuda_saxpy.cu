// Minimal CUDA smoke kernel: SAXPY (y = a*x + y).
// Used to verify the nvcc -> g++ link path end-to-end before any real
// arithmetic is ported. Only compiled when nvcc is detected by Makevars.

#include <cuda_runtime.h>

__global__ void saxpy_kernel(int n, float a, const float* x, float* y) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) y[i] = a * x[i] + y[i];
}

extern "C" int didgpu_cuda_saxpy(int n, float a, const float* x, float* y) {
  float *dx = nullptr, *dy = nullptr;
  size_t bytes = static_cast<size_t>(n) * sizeof(float);
  cudaError_t err;

  err = cudaMalloc(&dx, bytes);                                   if (err) goto done;
  err = cudaMalloc(&dy, bytes);                                   if (err) goto done;
  err = cudaMemcpy(dx, x, bytes, cudaMemcpyHostToDevice);          if (err) goto done;
  err = cudaMemcpy(dy, y, bytes, cudaMemcpyHostToDevice);          if (err) goto done;

  {
    const int block = 256;
    const int grid  = (n + block - 1) / block;
    saxpy_kernel<<<grid, block>>>(n, a, dx, dy);
  }
  err = cudaDeviceSynchronize();                                   if (err) goto done;
  err = cudaMemcpy(y, dy, bytes, cudaMemcpyDeviceToHost);

done:
  if (dx) cudaFree(dx);
  if (dy) cudaFree(dy);
  return static_cast<int>(err);
}
