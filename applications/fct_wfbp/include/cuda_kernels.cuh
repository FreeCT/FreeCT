#pragma once
#ifndef __CUDA_KERNELS_H__
#define __CUDA_KERNELS_H__

__constant__ struct CTGeometry d_cg;
__constant__ struct ReconConfig d_rp;
__constant__ struct GPUPrecompute d_gpu_precompute;

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

  // ORG Implementation
  //int out_idx = channel_idx + proj_idx*d_cg.num_detector_cols;
  //output[out_idx] = tex2D(tex_row_sheet, beta_idx + 0.5f, alpha_idx + 0.5);
  
  int channel_offset = (d_cg.num_detector_cols_padded_fft - d_cg.num_detector_cols)/2;
  int out_idx = (channel_idx+channel_offset)  + proj_idx*(d_cg.num_detector_cols_padded_fft);
  output[out_idx] = tex2D(tex_row_sheet, beta_idx + 0.5f, alpha_idx + 0.5);
}

__global__ void multiply_filter(cufftComplex * row_sheet_fourier_domain, cufftComplex * filter){
  int channel_idx = threadIdx.x + blockDim.x*blockIdx.x;
  int proj_idx    = threadIdx.y + blockDim.y*blockIdx.y;

  if (channel_idx >= (d_cg.num_detector_cols_padded_fft/2 + 1))
    return;

  //int idx = channel_idx + proj_idx * (d_cg.num_detector_cols/2 +1);
  int idx = channel_idx + proj_idx * (d_cg.num_detector_cols_padded_fft/2 +1);

  float a = filter[channel_idx].x;
  float b = filter[channel_idx].y;
  float c = row_sheet_fourier_domain[idx].x;
  float d = row_sheet_fourier_domain[idx].y;
  
  row_sheet_fourier_domain[idx].x = (1.0f/d_cg.num_detector_cols_padded_fft)*(a*c - b*d); 
  row_sheet_fourier_domain[idx].y = (1.0f/d_cg.num_detector_cols_padded_fft)*(a*d + b*c);
}

__global__ void reshape_rebin_into_final_array(float * d_final_projection_array, float * d_filtered_row_sheet, int row_idx){
  //int channel_idx = threadIdx.x + blockDim.x*blockIdx.x;
  //int proj_idx    = threadIdx.y + blockDim.y*blockIdx.y;
  //
  //int in_idx = channel_idx + proj_idx * d_cg.num_detector_cols;
  //int out_idx = channel_idx + row_idx * d_cg.num_detector_cols + proj_idx * d_cg.num_detector_cols * d_cg.num_detector_rows;
  //
  //d_final_projection_array[out_idx] = d_filtered_row_sheet[in_idx];

  int channel_idx = threadIdx.x + blockDim.x*blockIdx.x;
  int proj_idx    = threadIdx.y + blockDim.y*blockIdx.y;

  int channel_offset = (d_cg.num_detector_cols_padded_fft - d_cg.num_detector_cols)/2;
  int in_idx = (channel_idx + channel_offset) + (proj_idx * d_cg.num_detector_cols_padded_fft);
  int out_idx = channel_idx + row_idx * d_cg.num_detector_cols + proj_idx * d_cg.num_detector_cols * d_cg.num_detector_rows;

  d_final_projection_array[out_idx] = d_filtered_row_sheet[in_idx];
}

#define Q 0.6f
__device__ float W(float q){
  float wq;
  q = fabs(q);
  if (q<Q)
    wq = 1.0f;
  else if (Q<=q && q<1.0f){
    wq = cosf((0.5f*3.14159265f)*(q-Q)/(1.0f-Q));
    wq = wq*wq;
  }
  else
    wq = 0.0f;

  return wq;
}

__global__ void backproject_kernel_naive(float * d_projection_data,
                                         float * d_reconstruction_data,
                                         float * d_tube_angles,
                                         float * d_table_positions,
                                         float * d_slice_locations){

  // This kernel is kept here for illustrative purposes, however it not utilized.
  // It is more verbose, and thus easier to understand, however less efficient than
  // backproject_kernel.

  // Figure out which pixel we're operating on
  int x_idx = threadIdx.x + blockDim.x*blockIdx.x;
  int y_idx = threadIdx.y + blockDim.y*blockIdx.y;
  int slice_idx = threadIdx.z + blockDim.z*blockIdx.z;

  float x_loc = (x_idx - 0.5f*(d_rp.nx -1.0f))*d_rp.recon_fov/d_rp.nx + d_rp.x_origin; // Could be moved to a memory-operation (i.e. reduce FMA)
  float y_loc = (y_idx - 0.5f*(d_rp.ny -1.0f))*d_rp.recon_fov/d_rp.ny + d_rp.y_origin; // Could be moved to a memory-operation (i.e. reduce FMA)
  float z_loc = d_slice_locations[slice_idx]; // GLOBAL READ

  // Precalculate a "constants" structure to reduce the number of multiplications required in kernel
  // 0.5*acquisition_fov = half_acquisition_fov;
  // 0.5*projections_per_rotation = projections_per_half_turn;
  // floorf(d_cg.total_number_of_projections/projections_per_half_turn) = n_half_turns
  // pow(distance_source_to_isocenter,2.0) = distance_source_to_isocenter_squared
  // z_rot/(2*pi) = z_rot_over_2_pi
  // 1/distance_source_to_detector = reciprocal_distance_source_to_detector
  // tanf(0.5*theta_cone) = tanf_theta_cone
  // d_cg.distance_source_to_detector/(d_cg.distance_source_to_isocenter*d_cg.detector_pixel_size_col) =  pixel_scale 
  
  // If outside of the acquisition FOV, bail out
  if (x_loc*x_loc + y_loc*y_loc > (0.5f*d_cg.acquisition_field_of_view)*(0.5f*d_cg.acquisition_field_of_view))
    return;

  int projections_per_half_turn = (d_cg.projections_per_rotation/2);
  int n_half_turns = floorf(d_cg.total_number_of_projections/projections_per_half_turn);
  
  for (int theta_tilde_idx = 0; theta_tilde_idx < projections_per_half_turn; theta_tilde_idx++){
    float backprojected_value  = 0.0f;
    float normalization_factor = 0.0f;
  
    for (int k = 0; k < n_half_turns; k++){
  
      int curr_theta_idx = theta_tilde_idx + k*(projections_per_half_turn);
      float curr_theta = d_tube_angles[curr_theta_idx]; // GLOBAL READ
      float curr_table_position = d_table_positions[curr_theta_idx]; // GLOBAL READ
  
      float p_hat = x_loc*sinf(curr_theta) - y_loc*cosf(curr_theta);
      float l_hat = sqrtf(powf(d_cg.distance_source_to_isocenter,2.0f) - powf(p_hat,2.0f)) - (x_loc * cosf(curr_theta) - y_loc * sinf(curr_theta));
      float q_hat = (z_loc - curr_table_position + (d_cg.z_rot/(2.0f*3.14159265f))*asinf(p_hat/d_cg.distance_source_to_isocenter))/(l_hat*tanf(0.5f*d_cg.theta_cone));
  
      // Current projection has a ray that intersects the voxel
      // (Unfortunately I don't think there is a way to avoid the condition branch here)
      if (q_hat > 1.0f || q_hat < -1.0f)
        continue;
  
      // Compute the interpolation "metadata"
      float pixel_scale = d_cg.distance_source_to_detector/(d_cg.distance_source_to_isocenter*d_cg.detector_pixel_size_col);
  
      float p_idx = p_hat*pixel_scale + d_cg.detector_central_col;
      float q_idx = 0.5f*(q_hat + 1.0f)*(d_cg.num_detector_rows-1.0f);
      
      int q_floor = floorf(q_idx);
      int q_ceil = ceilf(q_idx);
  
      int p_floor = floorf(p_idx);
      int p_ceil = ceilf(p_idx);
  
      if (p_floor < 0 || q_floor < 0)
        continue;
  
      if (p_ceil >= d_cg.num_detector_cols-1 || q_ceil >= d_cg.num_detector_rows-1)
        continue;
  
      float w_q = q_idx - (float)q_floor;
      float w_p = p_idx - (float)p_floor;
  
      // Carry out the interpolation
      int lookup_idx_low  = p_floor +   (q_floor * d_cg.num_detector_cols)   + (curr_theta_idx * d_cg.num_detector_rows * d_cg.num_detector_cols);
      int lookup_idx_high = p_floor + ((q_floor+1) * d_cg.num_detector_cols) + (curr_theta_idx * d_cg.num_detector_rows * d_cg.num_detector_cols);
      float interpolated_value =
        (1.0f - w_q)*(1.0f - w_p) * d_projection_data[lookup_idx_low] + (1.0f - w_q)*(w_p) * d_projection_data[lookup_idx_low+1] +
        (w_q)*(1.0f - w_p) * d_projection_data[lookup_idx_high]             + (w_q)*(w_p)  * d_projection_data[lookup_idx_high+1];
  
  
      float WEIGHT = W(q_hat);
      backprojected_value += WEIGHT*interpolated_value;
      normalization_factor += WEIGHT;
    }
    
    // After looping over all half angles, add the result back to the slice data
    int pixel_idx = x_idx + (y_idx * d_rp.nx) + (slice_idx * d_rp.nx * d_rp.ny);
    d_reconstruction_data[pixel_idx] += (2.0f*3.1415926f/d_cg.projections_per_rotation)*(1.0f/normalization_factor) * backprojected_value;
    
  }
  
  
  // Done  
  //int out_idx = x_idx  + (y_idx * d_rp.nx) + (slice_idx * d_rp.nx * d_rp.ny);
  //d_reconstruction_data[out_idx] = z_loc;
  //d_reconstruction_data[out_idx] = 0.0f;
  
}

__global__ void backproject_kernel(float * d_projection_data,
                                   float * d_reconstruction_data,
                                   float * d_tube_angles,
                                   float * d_table_positions,
                                   float * d_slice_locations){

  // Figure out which pixel we're operating on
  int x_idx = threadIdx.x + blockDim.x*blockIdx.x;
  int y_idx = threadIdx.y + blockDim.y*blockIdx.y;
  int slice_idx = threadIdx.z + blockDim.z*blockIdx.z;

  // NOTE: This usage of shared memory neither helps, nor hurts
  // Unclear if we should bother keeping it around
  extern __shared__ float f[];
  float * x_data = f;
  float * y_data = f + d_rp.nx;
  
  x_data[x_idx] = (x_idx - 0.5f*(d_rp.nx -1.0f))*d_rp.recon_fov/d_rp.nx + d_rp.x_origin; // Could be moved to a memory-operation (i.e. reduce FMA)
  y_data[y_idx] = (y_idx - 0.5f*(d_rp.ny -1.0f))*d_rp.recon_fov/d_rp.ny + d_rp.y_origin; // Could be moved to a memory-operation (i.e. reduce FMA)
  float z_loc = d_slice_locations[slice_idx];
  
  // If outside of the acquisition FOV, bail out
  if (x_data[x_idx]*x_data[x_idx] + y_data[y_idx]*y_data[y_idx] > d_gpu_precompute.half_acquisition_fov_squared)
    return;

  for (int theta_tilde_idx = 0; theta_tilde_idx < d_gpu_precompute.projections_per_half_turn; theta_tilde_idx++){
    float backprojected_value  = 0.0f;
    float normalization_factor = 0.0f;
  
    for (int k = 0; k < d_gpu_precompute.n_half_turns; k++){
  
      int curr_theta_idx = theta_tilde_idx + k*(d_gpu_precompute.projections_per_half_turn);
      float curr_table_position = d_table_positions[curr_theta_idx]; // GLOBAL READ

      if (fabsf(curr_table_position - z_loc) > d_cg.z_rot)
        continue;

      float curr_theta = d_tube_angles[curr_theta_idx]; // GLOBAL READ
      
      float p_hat = x_data[x_idx]*sinf(curr_theta) - y_data[y_idx]*cosf(curr_theta);
      float l_hat = sqrtf(d_gpu_precompute.distance_source_to_isocenter_squared - powf(p_hat,2.0f)) - x_data[x_idx] * cosf(curr_theta) - y_data[y_idx] * sinf(curr_theta);
      float q_hat = (z_loc - curr_table_position + (d_gpu_precompute.z_rot_over_2_pi)*asinf(d_gpu_precompute.recip_distance_source_to_isocenter*p_hat))/(l_hat*d_gpu_precompute.tanf_theta_cone);
      
      // Current projection has a ray that intersects the voxel
      // (Unfortunately I don't think there is a way to avoid the condition branch here)
      if (q_hat > 1.0f || q_hat < -1.0f)
        continue;
  
      // Compute the interpolation "metadata"      
      float p_idx = p_hat*d_gpu_precompute.pixel_scale + d_cg.detector_central_col;
      float q_idx = 0.5*(q_hat + 1.0f)*(d_cg.num_detector_rows-1.0f);

      //p_idx = (d_cg.num_detector_cols+1.5) - p_idx;
      p_idx = (d_cg.num_detector_cols+1.75) - p_idx;
      
      int q_floor = floorf(q_idx);
      int q_ceil = ceilf(q_idx);
  
      int p_floor = floorf(p_idx);
      int p_ceil = ceilf(p_idx);
  
      if (p_floor < 0 || q_floor < 0)
        continue;
  
      if (p_ceil >= d_cg.num_detector_cols-1 || q_ceil >= d_cg.num_detector_rows-1)
        continue;
  
      float w_q = q_idx - (float)q_floor;
      float w_p = p_idx - (float)p_floor;
  
      // Carry out the interpolation
      int lookup_idx_low  = p_floor +   (q_floor * d_cg.num_detector_cols)   + (curr_theta_idx * d_cg.num_detector_rows * d_cg.num_detector_cols);
      int lookup_idx_high = p_floor + ((q_floor+1) * d_cg.num_detector_cols) + (curr_theta_idx * d_cg.num_detector_rows * d_cg.num_detector_cols);
      float interpolated_value =
        (1.0f - w_q)*(1.0f - w_p) * d_projection_data[lookup_idx_low] + (1.0f - w_q)*(w_p) * d_projection_data[lookup_idx_low+1] +
        (w_q)*(1.0f - w_p) * d_projection_data[lookup_idx_high]             + (w_q)*(w_p)  * d_projection_data[lookup_idx_high+1];
  
  
      float WEIGHT = W(q_hat);
      backprojected_value += WEIGHT*interpolated_value;
      normalization_factor += WEIGHT; 
    }
    
    // After looping over all half angles, add the result back to the slice data
    int pixel_idx = x_idx + (y_idx * d_rp.nx) + (slice_idx * d_rp.nx * d_rp.ny);
    d_reconstruction_data[pixel_idx] += d_gpu_precompute.two_pi_over_projections_per_rotation * (1.0f/normalization_factor) * backprojected_value;
    
  }
  
}

__global__ void thicken_slices(float * d_recon_native_thickness,
                               float * d_recon_requested_thickness,
                               float * d_slice_locations_native,
                               float * d_slice_locations_requested){

  // Index over the ouput array
  int x_idx = threadIdx.x + blockDim.x*blockIdx.x;
  int y_idx = threadIdx.y + blockDim.y*blockIdx.y;
  int slice_idx = threadIdx.z + blockDim.z*blockIdx.z;

  float curr_slice_location = d_slice_locations_requested[slice_idx];

  float pixel_val = 0.0f;
  float weight_sum = 0.0f;
  for (int i=0;i<d_gpu_precompute.n_slices_native;i++){
  
    float weight = fmaxf(0.0f,1.0f - fabsf(d_slice_locations_native[i] - curr_slice_location)/d_rp.slice_thickness);
    
    if (weight==0.0f)
      continue;
    
    weight_sum += weight;
    
    int in_idx = x_idx + (y_idx * d_rp.nx) + (i * d_rp.nx * d_rp.ny);
    pixel_val += weight * d_recon_native_thickness[in_idx];
  }

  int out_idx = x_idx + (y_idx * d_rp.nx) + (slice_idx * d_rp.nx * d_rp.ny);
  d_recon_requested_thickness[out_idx] = pixel_val/weight_sum;
  //d_recon_requested_thickness[out_idx] = pixel_val/d_gpu_precompute.n_slices_native;
  
}


#endif // __CUDA_KERNELS_H__