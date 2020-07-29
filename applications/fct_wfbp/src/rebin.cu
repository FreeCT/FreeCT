#include <rebin.h>
#include <iostream>
#include <util.h>

#include <gpu_error_check.h>
#include <cufft.h>

#include <chrono>

// @JOHN: Don't sweat seconds: focus on clean code and algorithmic improvements

// GPU CODE block
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
  //output[out_idx] = d_cg.detector_pixel_size_col; //tex2D(tex_row_sheet, channel_idx + 0.5f, proj_idx + 0.5);

}

__global__ void multiply_filter(cufftComplex * row_sheet_fourier_domain, float * filter){

  int channel_idx = threadIdx.x + blockDim.x*blockIdx.x;
  int proj_idx    = threadIdx.y + blockDim.y*blockIdx.y;

  int idx = channel_idx + proj_idx * d_cg.num_detector_cols;
  row_sheet_fourier_domain[idx] = filter[channel_idx] * row_sheet_fourier_domain[channel_idx];

}

inline void configure_texture(){
  tex_row_sheet.addressMode[0] = cudaAddressModeClamp;
  tex_row_sheet.addressMode[1] = cudaAddressModeClamp;
  tex_row_sheet.addressMode[2] = cudaAddressModeClamp;
  tex_row_sheet.filterMode     = cudaFilterModeLinear;
  tex_row_sheet.normalized     = false;
}

// HOST CODE
void generate_filter(float * f_array,CTGeometry cg, float c = 1.0f, float a = 1.0f);
  
void rebin(std::shared_ptr<float> output, std::shared_ptr<float> input, CTGeometry cg, ReconConfig rp){

  // Reshape the raw data to be "row sheets"
  // Additionally, flip the channel direction and row direction
  // since Chen et al 2015 defines the geometry to be the opposite
  // of how Stierstorfer et al 2004 defines it.
  Timer t;
  t.tic();
  std::cout << "Reshaping raw data array..." << std::endl;
  std::shared_ptr<float> raw_reshaped_ptr(new float[cg.num_detector_cols*cg.num_detector_rows*cg.total_number_of_projections]);
  float * raw = input.get();
  float * raw_reshaped = raw_reshaped_ptr.get();
  for (int i=0; i<cg.total_number_of_projections;i++){
    for (int j=0; j<cg.num_detector_rows;j++){
      for (int k=0; k<cg.num_detector_cols;k++){
        int input_idx = k + j*cg.num_detector_cols + i*cg.num_detector_cols*cg.num_detector_rows;
        int output_idx = (cg.num_detector_cols - 1 - k) + i*cg.num_detector_cols + (cg.num_detector_rows - 1 - j)*cg.num_detector_cols*cg.total_number_of_projections;
        raw_reshaped[output_idx] = raw[input_idx];
      }
    }
  }
  std::cout << "Done!" << std::endl;
  t.toc();

  // Allocate our GPU arrays
  cudaError_t gpu_status;
  
  cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc<float>();
  cudaArray * d_row_sheet_raw;
  gpu_status = cudaMallocArray(&d_row_sheet_raw, &channelDesc, cg.num_detector_cols, cg.total_number_of_projections);
  gpuErrChk(gpu_status);
  
  float * d_row_sheet_rebin;
  gpu_status = cudaMalloc(&d_row_sheet_rebin,cg.num_detector_cols*cg.total_number_of_projections*sizeof(float));
  gpuErrChk(gpu_status);
  
  gpu_status = cudaMemcpyToSymbol(d_cg,&cg,sizeof(struct CTGeometry),0,cudaMemcpyHostToDevice);
  gpuErrChk(gpu_status);
  
  gpu_status = cudaMemcpyToSymbol(d_rp,&rp,sizeof(struct ReconConfig),0,cudaMemcpyHostToDevice);
  gpuErrChk(gpu_status);

  configure_texture();

  // Create our filter and send to GPU
  float * h_filter = new float[]
  
  
  // Main rebin/filter loop
  float * rebinned_data = output.get();
  
  GPUTimer gt;
  for (int i=0; i<cg.num_detector_rows; i++){
    gt.tic();    
    size_t offset = i*cg.num_detector_cols*cg.total_number_of_projections;
    size_t sheet_size_bytes =  cg.num_detector_cols*cg.total_number_of_projections * sizeof(float);
    gpu_status = cudaMemcpyToArray(d_row_sheet_raw, 0, 0, &raw_reshaped[offset], sheet_size_bytes, cudaMemcpyHostToDevice);
    gpuErrChk(gpu_status);

    gpu_status = cudaBindTextureToArray(tex_row_sheet,d_row_sheet_raw,channelDesc);
    gpuErrChk(gpu_status);

    dim3 rebin_threads(8,8);
    dim3 rebin_blocks(cg.num_detector_cols/rebin_threads.x,cg.total_number_of_projections/rebin_threads.y);
    rebin_kernel<<<rebin_blocks,rebin_threads>>>(d_row_sheet_rebin);
    gpuErrChk(cudaPeekAtLastError());

    // Configure FFT plan and run filtering
    cufftHandle plan;
    cufftComplex * d_sheet_data_fourier_domain;

    int RANK = 1;
    int NX = cg.num_detector_cols;
    int BATCH = cg.total_number_of_projections;
    
    cudaMalloc(&d_sheet_data_fourier_domain,cg.num_detector_cols*cg.total_number_of_projections*sizeof(cufftComplex));

    cufftPlanMany(&plan,RANK,&NX,
                  NULL,1,0,
                  NULL,1,0,
                  CUFFT_R2C,BATCH);
    cufftExecR2C(plan,d_row_sheet_rebin,d_sheet_data_fourier_domain);

    dim3 filter_threads(cg.num_detector_cols,1);
    dim3 filter_blocks(1,cg.total_number_of_projections/rebin_threads.y);
    multiply_filter<<<filter_blocks,filter_threads>>>(d_sheet_data_fourier_domain,);
    
    // Copy data back from GPU
    cudaMemcpy(&rebinned_data[offset],d_row_sheet_rebin,sheet_size_bytes,cudaMemcpyDeviceToHost);


      
    gt.toc();
    
  }
  
  cudaFreeArray(d_row_sheet_raw);
  cudaFree(d_row_sheet_rebin);
}

void generate_filter(float * f_array, CTGeometry cg, float c, float a){
  // Create a spatial domain ramp filter.  Eventually we'll expost c and a
  // so users can customize filter response for smoother/sharper reconstructions
  
  //float * h_filter=(float*)calloc(2*cg.n_channels_oversampled,sizeof(float));
  //float ds = mr->cg.r_f*sin(mr->cg.fan_angle_increment/2.0f); // This is at isocenter.  Is that correct?

  float pi_f = 3.141592653589f;
  float ds = cg.detector_pixel_size_col;
  
  //float ds = mr->cg.src_to_det*sin(mr->cg.fan_angle_increment/2.0f); // I think it should be at the detector

  auto r = [](float t)-> float{
             float v = sin(t)/t + (cos(t)-1.0f)/(t*t);
             if (t==0)
               v=0.5;
             return v;
           };
  
  int test = (int)cg.num_detector_cols;
  
  for (int i = -test;i < test;i++){
    
    f_array[i+cg.num_detector_cols] = (c*c/(2.0f*ds)) * (a*r(c*pi_f*i) +
                                                         (((1.0f-a)/2.0f)*r(pi_f*c*i + pi_f)) +
                                                         (((1.0f -a)/2.0f)*r(pi_f*c*i-pi_f)));
    
  }
  
}
