#pragma once
#ifndef __CUDA_KERNELS_H__
#define __CUDA_KERNELS_H__

__constant__ struct CTGeometry d_cg;
__constant__ struct ReconConfig d_rp;
texture<float,cudaTextureType2D,cudaReadModeElementType> tex_row_sheet;

__global__ void rebin_kernel(float * output){
  int channel_idx = threadIdx.x + blockDim.x*blockIdx.x;
  int proj_idx    = threadIdx.y + blockDim.y*blockIdx.y;

  float delta_theta = 2.0f*3.14159265f/(float)d_cg.projections_per_rotation;
  
  float p = ((float)channel_idx - d_cg.detector_central_col)*(float)d_cg.detector_pixel_size_col;
  float beta = asin(p/d_cg.distance_source_to_detector);
  float theta = (float)proj_idx * delta_theta;
  float alpha = theta - beta;
  
  float beta_idx  = (beta*d_cg.distance_source_to_detector/d_cg.detector_pixel_size_col) + d_cg.detector_central_col;
  float alpha_idx = alpha/delta_theta;
  
  int out_idx = channel_idx + proj_idx*d_cg.num_detector_cols;
  output[out_idx] = tex2D(tex_row_sheet, beta_idx + 0.5f, alpha_idx + 0.5);
}

__global__ void multiply_filter(cufftComplex * row_sheet_fourier_domain, cufftComplex * filter){
  int channel_idx = threadIdx.x + blockDim.x*blockIdx.x;
  int proj_idx    = threadIdx.y + blockDim.y*blockIdx.y;

  int idx = channel_idx + proj_idx * (d_cg.num_detector_cols/2 +1);

  float a = filter[channel_idx].x;
  float b = filter[channel_idx].y;
  float c = row_sheet_fourier_domain[idx].x;
  float d = row_sheet_fourier_domain[idx].y;
  
  row_sheet_fourier_domain[idx].x = a*c - b*d;
  row_sheet_fourier_domain[idx].y = a*d + b*c;
}

__global__ void reshape_rebin_into_final_array(float * d_final_projection_array, float * d_filtered_row_sheet, int row_idx){
  int channel_idx = threadIdx.x + blockDim.x*blockIdx.x;
  int proj_idx    = threadIdx.y + blockDim.y*blockIdx.y;

  int in_idx = channel_idx + proj_idx * d_cg.num_detector_cols;
  int out_idx = channel_idx + row_idx * d_cg.num_detector_cols + proj_idx * d_cg.num_detector_cols * d_cg.num_detector_rows;

  d_final_projection_array[out_idx] = d_filtered_row_sheet[in_idx];
}

#endif // __CUDA_KERNELS_H__