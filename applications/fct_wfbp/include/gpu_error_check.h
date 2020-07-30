#pragma once

#include <cuda_runtime.h>
#include <cufft.h>

class GPUTimer{
public:
  GPUTimer(){
    cudaEventCreate(&m_start);
    cudaEventCreate(&m_stop);
  };
  ~GPUTimer(){
    cudaEventDestroy(m_start);
    cudaEventDestroy(m_stop);
  };

  void tic(){
    cudaEventRecord(m_start);
  }

  void toc(){
    cudaEventRecord(m_stop);
    cudaEventSynchronize(m_stop);
    
    float milliseconds = 0.0f;
    cudaEventElapsedTime(&milliseconds,m_start,m_stop);
    std::cout << "GPU timer: " << milliseconds << "ms" << std::endl;
  }

private:
  cudaEvent_t m_start;
  cudaEvent_t m_stop;

};

// CUDA runtime error checking

#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess) 
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
    }
}


// cufft error checking
static const char * cufftGetErrorString(cufftResult error);

#define cufftErrChk(ans){cufftAssert((ans),__FILE__,__LINE__);}
inline void cufftAssert(cufftResult code, const char *file, int line, bool abort=true){
  if (code != CUFFT_SUCCESS) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cufftGetErrorString(code), file, line);
      if (abort) exit(code);
  }
}

static const char *cufftGetErrorString(cufftResult error)
{
  switch (error)
    {
    case CUFFT_SUCCESS:
      return "CUFFT_SUCCESS";

    case CUFFT_INVALID_PLAN:
      return "CUFFT_INVALID_PLAN";

    case CUFFT_ALLOC_FAILED:
      return "CUFFT_ALLOC_FAILED";

    case CUFFT_INVALID_TYPE:
      return "CUFFT_INVALID_TYPE";

    case CUFFT_INVALID_VALUE:
      return "CUFFT_INVALID_VALUE";

    case CUFFT_INTERNAL_ERROR:
      return "CUFFT_INTERNAL_ERROR";

    case CUFFT_EXEC_FAILED:
      return "CUFFT_EXEC_FAILED";

    case CUFFT_SETUP_FAILED:
      return "CUFFT_SETUP_FAILED";

    case CUFFT_INVALID_SIZE:
      return "CUFFT_INVALID_SIZE";

    case CUFFT_UNALIGNED_DATA:
      return "CUFFT_UNALIGNED_DATA";
    }

  return "<unknown>";
}