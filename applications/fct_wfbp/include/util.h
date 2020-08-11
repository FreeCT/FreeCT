#pragma once

#include <chrono>
#include <iostream>

class Timer{
public:
  Timer(){};
  ~Timer(){};

  void tic(){
    m_start = std::chrono::high_resolution_clock::now();
  }
  void toc(){
    m_end = std::chrono::high_resolution_clock::now();
    std::cout << "Required "<< (float)std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start).count()/1000.0f << "s" << std::endl;    
  }
private:
  std::chrono::high_resolution_clock::time_point m_start;
  std::chrono::high_resolution_clock::time_point m_end;
};

bool validate_selected_device(int device_idx);

void debug_save_cuda_array(float * device_array, size_t num_bytes, std::string filepath);