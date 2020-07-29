#pragma once

#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess) 
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
    }
}

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