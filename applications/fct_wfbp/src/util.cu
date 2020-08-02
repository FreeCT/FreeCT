#include <util.h>
#include <iostream>

bool validate_selected_device(int cuda_device){

  // Ensure the requested device even makes sense
  int device_count=0;
  cudaGetDeviceCount(&device_count);
  if (device_count < 1){
    std::cout << "ERROR: No CUDA devices detected." << std::endl;
    return false;
  }
  
  if (cuda_device > device_count - 1){
    std::cout << "ERROR: Invalid CUDA device requested!" << std::endl;
    std::cout << "       Requested value (" << cuda_device << ") must be less than the number of devices (" << device_count  << ")" << std::endl;
    return false;
  }

  // Attempt to set the cuda device
  cudaSetDevice(cuda_device);

  // Double check that the device was set correctly
  int final_idx;
  cudaGetDevice(&final_idx);
  if (final_idx!=cuda_device){
    std::cout << "ERROR: Attempted to set the CUDA device however something has gone wrong." << std::endl;
    std::cout << "       Set the device to \"" << cuda_device << "\", however CUDA is reporting device as \"" << final_idx << "\n" << std::endl;
    return false;
  }
  
  return true;
}