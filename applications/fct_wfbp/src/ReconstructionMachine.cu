#include <ReconstructionMachine.h>

#include <fstream>

#include <boost/filesystem.hpp>

#include <gpu_error_check.h>
#include <parse_config.h>
#include <util.h>

#include <cuda_kernels.cuh>

#include <cufft.h>

namespace{
  cufftComplex * generate_filter(CTGeometry cg, float c = 1.0f, float a = 1.0f);
}

namespace fct{

  void ReconstructionMachine::LoadReconConfiguration(std::string filepath){
    std::cout << "ReconstructionMachine: Loading reconstruction configuration..." << std::endl;
    if (!boost::filesystem::exists(filepath)){
      std::cout << "ERROR: Recon configuration filepath does not exist! (" << filepath << ")" << std::endl;
      exit(1);
    }

    parse_config(filepath,m_rp);
  }
  
  void ReconstructionMachine::LoadRawData(){
    std::cout << "ReconstructionMachine: Loading raw data..." << std::endl;
    m_org_data_set = std::make_shared<fct::DicomDataSet>();
    m_org_data_set->setPath(m_rp.raw_data_dir);
    m_org_data_set->initialize();
    m_org_data_set->readAll();
  }
  
  void ReconstructionMachine::ConfigureCTGeometry(){
    std::cout << "ReconstructionMachine: Configuring CT geometry... " << std::endl;
    
    // Physical geometry of the scanner (cannot change from scan to scan)
    m_cg.total_number_of_projections  = m_org_data_set->getTotalNumProjections();
    m_cg.projections_per_rotation     = m_org_data_set->getProjectionsPerRotation();
    m_cg.detector_pixel_size_col     = m_org_data_set->getDetectorTransverseSpacing();
    m_cg.detector_pixel_size_row     = m_org_data_set->getDetectorAxialSpacing();
    m_cg.num_detector_cols            = m_org_data_set->getDetectorChannels();
    m_cg.num_detector_rows            = m_org_data_set->getDetectorRows();
    m_cg.detector_central_col         = m_org_data_set->getDetectorCentralChannel();
    m_cg.detector_central_row         = m_org_data_set->getDetectorCentralRow();
    m_cg.distance_source_to_isocenter = m_org_data_set->getDistSourceToIsocenter(); 
    m_cg.distance_source_to_detector  = m_org_data_set->getDistSourceToDetector();
    
    m_cg.collimated_slice_width       = m_org_data_set->getDistSourceToIsocenter()*(m_org_data_set->getDetectorAxialSpacing()/m_org_data_set->getDistSourceToDetector());
    
    m_cg.z_rot = fabs(m_org_data_set->getTablePosition(m_cg.projections_per_rotation) - m_org_data_set->getTablePosition(0));
  
    float detector_cone_offset = ((float)(m_cg.num_detector_rows - 1))/2.0f; // May not be 100% accurate if central detector is not necessarily in the middle
    m_cg.theta_cone=2.0f*atan(detector_cone_offset * m_cg.collimated_slice_width/m_cg.distance_source_to_isocenter);

    m_cg.acquisition_field_of_view = 2.0f * m_cg.distance_source_to_isocenter*sin((float(m_cg.num_detector_cols-1.0f)/2.0f) * m_org_data_set->getDetectorTransverseSpacing() * (1.0f/m_cg.distance_source_to_detector));

    m_gpu_precompute.InitFromCTGeometry(m_cg);
  }
  
  void ReconstructionMachine::RunReconstruction(){
    std::cout << "ReconstructionMachine: Running the reconstruction (INCOMPLETE!)..." << std::endl;
    AllocateKeyArrays();
    RebinAndFilter();
    Backproject();
    ApplyFinalSliceThicknessing();
  }

  void ReconstructionMachine::AllocateKeyArrays(){

    /* HOST ARRAYS */
    // Set up the raw data array suitable for CUDA
    m_raw_data.reset(new float[m_cg.num_detector_rows*m_cg.num_detector_cols*m_cg.total_number_of_projections]);

    for (size_t projection_idx=0; projection_idx < m_cg.total_number_of_projections;projection_idx++){
      size_t data_offset = m_cg.num_detector_cols * m_cg.num_detector_rows * projection_idx;
      m_org_data_set->copyProjection(projection_idx,&(m_raw_data.get()[data_offset]));
    }

    // Tube angles and table positions
    m_tube_angles.resize(m_cg.total_number_of_projections);
    m_table_positions.resize(m_cg.total_number_of_projections);
    
    for (size_t projection_idx=0; projection_idx < m_cg.total_number_of_projections;projection_idx++){
      m_table_positions[projection_idx] = m_org_data_set->getTablePosition(projection_idx);
      m_tube_angles[projection_idx] = m_org_data_set->getTubeAngle(projection_idx);
    }
    
    // reclaim some memory
    m_org_data_set.reset();

    // Compute the locations of the slices (both in the native
    // detector slice thickness as well as at the requested slice
    // pitch.
    if (m_rp.start_pos!=m_rp.end_pos){
      float slice_delta_native = m_cg.collimated_slice_width * fabs(m_rp.end_pos - m_rp.start_pos)/(m_rp.end_pos - m_rp.start_pos);
      float slice_delta = m_rp.slice_pitch * fabs(m_rp.end_pos - m_rp.start_pos)/(m_rp.end_pos - m_rp.start_pos);
      
      if (slice_delta > 0){
        for (float loc = m_rp.start_pos; loc<m_rp.end_pos; loc+=slice_delta_native)
          m_slice_locations_collimated_slice_width.push_back(loc);

        for (float loc = m_rp.start_pos; loc<m_rp.end_pos; loc+=slice_delta)
          m_slice_locations_requested_slice_width.push_back(loc);
      }
      else{
        for (float loc = m_rp.start_pos; loc>m_rp.end_pos; loc+=slice_delta_native)
          m_slice_locations_collimated_slice_width.push_back(loc);

        for (float loc = m_rp.start_pos; loc>m_rp.end_pos; loc+=slice_delta)
          m_slice_locations_requested_slice_width.push_back(loc);        
      }
    }
    else{
      m_slice_locations_collimated_slice_width.push_back(m_rp.start_pos);
      m_slice_locations_requested_slice_width.push_back(m_rp.start_pos);      
    }
  }
  
  void ReconstructionMachine::RebinAndFilter(){
    // Reshape the raw data to be "row sheets"
    // Additionally, flip the channel direction and row direction
    // since Chen et al 2015 defines the geometry to be the opposite
    // of how Stierstorfer et al 2004 defines it.
    // (May be worth it to eventually do this on the GPU...)
    Timer t;
    t.tic();
    std::cout << "Reshaping raw data array..." << std::endl;
    std::shared_ptr<float> raw_reshaped_ptr(new float[m_cg.num_detector_cols*m_cg.num_detector_rows*m_cg.total_number_of_projections]);
    float * raw = m_raw_data.get();
    float * raw_reshaped = raw_reshaped_ptr.get();
    for (int i=0; i<m_cg.total_number_of_projections;i++){
      for (int j=0; j<m_cg.num_detector_rows;j++){
        for (int k=0; k<m_cg.num_detector_cols;k++){
          int input_idx = k + j*m_cg.num_detector_cols + i*m_cg.num_detector_cols*m_cg.num_detector_rows;
          int output_idx = (m_cg.num_detector_cols - 1 - k) + i*m_cg.num_detector_cols + (m_cg.num_detector_rows - 1 - j)*m_cg.num_detector_cols*m_cg.total_number_of_projections;
          //int output_idx = (m_cg.num_detector_cols - 1 - k) + i*m_cg.num_detector_cols + (m_cg.num_detector_rows - 1 - j)*m_cg.num_detector_cols*m_cg.total_number_of_projections;
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
    gpu_status = cudaMallocArray(&d_row_sheet_raw, &channelDesc, m_cg.num_detector_cols, m_cg.total_number_of_projections);
    gpuErrChk(gpu_status);

    tex_row_sheet.addressMode[0] = cudaAddressModeClamp;
    tex_row_sheet.addressMode[1] = cudaAddressModeClamp;
    tex_row_sheet.addressMode[2] = cudaAddressModeClamp;
    tex_row_sheet.filterMode     = cudaFilterModeLinear;
    tex_row_sheet.normalized     = false;
  
    float * d_row_sheet_rebin;
    gpu_status = cudaMalloc(&d_row_sheet_rebin,m_cg.num_detector_cols*m_cg.total_number_of_projections*sizeof(float));
    gpuErrChk(gpu_status);
  
    gpu_status = cudaMemcpyToSymbol(d_cg,&m_cg,sizeof(struct CTGeometry),0,cudaMemcpyHostToDevice);
    gpuErrChk(gpu_status);
  
    gpu_status = cudaMemcpyToSymbol(d_rp,&m_rp,sizeof(struct ReconConfig),0,cudaMemcpyHostToDevice);
    gpuErrChk(gpu_status);

    gpu_status = cudaMemcpyToSymbol(d_gpu_precompute,&m_gpu_precompute,sizeof(struct GPUPrecompute),0,cudaMemcpyHostToDevice);
    gpuErrChk(gpu_status);

    // Allocate the array we'll reshape our final, rebinned, filtered data
    // into and pass to backprojection
    gpu_status  = cudaMalloc(&m_d_filtered_projection_data,m_cg.num_detector_cols*m_cg.num_detector_rows*m_cg.total_number_of_projections*sizeof(float));
    gpuErrChk(gpu_status);
    
    // Create our filter (on device, ready to be multiplied against our data)
    cufftComplex * d_filter = generate_filter(m_cg);

    // Create our FFT plans
    cufftResult cufft_status;
    cufftHandle plan_forward;
    cufftHandle plan_reverse;
  
    int RANK  = 1;
    int NX    = m_cg.num_detector_cols;
    int BATCH = m_cg.total_number_of_projections;

    cufft_status = cufftPlanMany(&plan_forward,RANK,&NX,
                                 NULL,1,0,
                                 NULL,1,0,
                                 CUFFT_R2C,BATCH);
    cufftErrChk(cufft_status);

    cufft_status = cufftPlanMany(&plan_reverse,RANK,&NX,
                                 NULL,1,0,
                                 NULL,1,0,
                                 CUFFT_C2R,BATCH);
    cufftErrChk(cufft_status);

    // Allocate FFT result data
    int fft_output_size = NX/2+1;
    cufftComplex * d_sheet_data_fourier_domain;
    gpu_status = cudaMalloc(&d_sheet_data_fourier_domain,fft_output_size*m_cg.total_number_of_projections*sizeof(cufftComplex));
    gpuErrChk(gpu_status);

    // Main rebin/filter loop
    std::cout << "Rebinning and filtering data..." << std::endl;

    //std::shared_ptr<float> output(new float[m_cg.num_detector_cols*m_cg.num_detector_rows*m_cg.total_number_of_projections]);
    //float * rebinned_data = output.get();
    
    GPUTimer gt;
    gt.tic();
    for (int i=0; i<m_cg.num_detector_rows; i++){
      
      // Copy raw projection data to texture memory
      size_t offset = i*m_cg.num_detector_cols*m_cg.total_number_of_projections;
      size_t sheet_size_bytes =  m_cg.num_detector_cols*m_cg.total_number_of_projections * sizeof(float);
      gpu_status = cudaMemcpyToArray(d_row_sheet_raw, 0, 0, &raw_reshaped[offset], sheet_size_bytes, cudaMemcpyHostToDevice);
      gpuErrChk(gpu_status);

      gpu_status = cudaBindTextureToArray(tex_row_sheet,d_row_sheet_raw,channelDesc);
      gpuErrChk(gpu_status);

      // Run rebinning kernel
      dim3 rebin_threads(8,8);
      dim3 rebin_blocks(m_cg.num_detector_cols/rebin_threads.x,m_cg.total_number_of_projections/rebin_threads.y);
      rebin_kernel<<<rebin_blocks,rebin_threads>>>(d_row_sheet_rebin);
      gpuErrChk(cudaPeekAtLastError());

      // Filter the rebinned data
      cufft_status = cufftExecR2C(plan_forward,(cufftReal*)d_row_sheet_rebin,d_sheet_data_fourier_domain);
      cufftErrChk(cufft_status);
      cudaDeviceSynchronize();

      dim3 filter_threads(fft_output_size,1);
      dim3 filter_blocks(1,m_cg.total_number_of_projections/filter_threads.y);
      multiply_filter<<<filter_blocks,filter_threads>>>(d_sheet_data_fourier_domain,d_filter);
      gpuErrChk(cudaPeekAtLastError());
      cudaDeviceSynchronize();

      cufft_status = cufftExecC2R(plan_reverse,d_sheet_data_fourier_domain,(cufftReal*)d_row_sheet_rebin);
      cufftErrChk(cufft_status);
      cudaDeviceSynchronize();
    
      // Copy data back from GPU
      //gpu_status = cudaMemcpy(&rebinned_data[offset],d_row_sheet_rebin,sheet_size_bytes,cudaMemcpyDeviceToHost);
      //gpuErrChk(gpu_status);

      // Reshape the data into the final projection array for backprojections
      dim3 reshape_threads(m_cg.num_detector_cols,1);
      dim3 reshape_blocks(1,m_cg.total_number_of_projections/filter_threads.y);
      reshape_rebin_into_final_array<<<reshape_blocks,reshape_threads>>>(m_d_filtered_projection_data,d_row_sheet_rebin,i);
      gpuErrChk(cudaPeekAtLastError());
    }
    gt.toc();

    cufftDestroy(plan_forward);
    cufftDestroy(plan_reverse);
  
    cudaFree(d_sheet_data_fourier_domain);
    cudaFree(d_filter);
    cudaFreeArray(d_row_sheet_raw);
    cudaFree(d_row_sheet_rebin);

    //std::ofstream fid("/home/john/Desktop/rebin_debug.dat",std::ios::binary);
    //fid.write((char*)rebinned_data,m_cg.num_detector_cols*m_cg.num_detector_rows*m_cg.total_number_of_projections*sizeof(float));
      
  }
  
  void ReconstructionMachine::Backproject(){

    std::cout << "Running backprojection for " << m_slice_locations_collimated_slice_width.size() << " slices ..." << std::endl;

    cudaError_t cuda_status;

    // Copy table positions and tube angles onto the gpu
    float * d_tube_angles;
    cuda_status = cudaMalloc(&d_tube_angles,m_tube_angles.size()*sizeof(float));
    gpuErrChk(cuda_status);
    cuda_status = cudaMemcpy(d_tube_angles,&m_tube_angles[0],m_tube_angles.size()*sizeof(float),cudaMemcpyHostToDevice);
    gpuErrChk(cuda_status);
    
    float * d_table_positions;
    cuda_status = cudaMalloc(&d_table_positions,m_table_positions.size()*sizeof(float));
    gpuErrChk(cuda_status);
    cuda_status = cudaMemcpy(d_table_positions,&m_table_positions[0],m_table_positions.size()*sizeof(float),cudaMemcpyHostToDevice);
    gpuErrChk(cuda_status);

    float * d_slice_locations_collimated_slice_width;
    cuda_status = cudaMalloc(&d_slice_locations_collimated_slice_width,m_slice_locations_collimated_slice_width.size()*sizeof(float));
    gpuErrChk(cuda_status);
    cuda_status = cudaMemcpy(d_slice_locations_collimated_slice_width,
                             &m_slice_locations_collimated_slice_width[0],
                             m_slice_locations_collimated_slice_width.size()*sizeof(float),
                             cudaMemcpyHostToDevice);
    gpuErrChk(cuda_status);
      
    // Allocate GPU memory for the reconstructed slices (at native slice thickness)
    int m_num_slices_native = m_slice_locations_collimated_slice_width.size();
    cuda_status = cudaMalloc(&m_d_reconstruction_collimated_slice_width, m_rp.nx * m_rp.ny * m_num_slices_native * sizeof(float));
    gpuErrChk(cuda_status);
    cuda_status = cudaMemset(m_d_reconstruction_collimated_slice_width,0,m_rp.nx * m_rp.ny * m_num_slices_native * sizeof(float));
    gpuErrChk(cuda_status);
    
    // All of our projection data is already on GPU and stored in m_d_filtered_projection_data from rebinning/filtering
    GPUTimer gt;
    gt.tic();
    dim3 backproject_threads(32,32,1);
    dim3 backproject_blocks(m_rp.nx / backproject_threads.x,
                            m_rp.ny / backproject_threads.y,
                            m_num_slices_native / backproject_threads.z);

    size_t shared_mem_size = m_rp.nx*sizeof(float) + m_rp.ny*sizeof(float);
    
    backproject_kernel<<<backproject_blocks,backproject_threads,shared_mem_size>>>(m_d_filtered_projection_data,
                                                                   m_d_reconstruction_collimated_slice_width,
                                                                   d_tube_angles,
                                                                   d_table_positions,
                                                                   d_slice_locations_collimated_slice_width);
    cudaDeviceSynchronize();
    gpuErrChk(cudaPeekAtLastError());
    gt.toc();
    
    // Debug
    if (true){
      m_reconstructed_data.reset(new float[m_rp.nx * m_rp.ny * m_num_slices_native]);

      if (m_reconstructed_data.get()==NULL){
        std::cout << "COULD NOT ALLOCATE RECONSTRUCTION DATA MEMORY" << std::endl;
        exit(1);
      }

      cuda_status = cudaMemcpy((void*)(m_reconstructed_data.get()),
                               (void*)m_d_reconstruction_collimated_slice_width,
                               m_rp.nx * m_rp.ny * m_num_slices_native*sizeof(float),cudaMemcpyDeviceToHost);
      gpuErrChk(cuda_status);
    
      std::ofstream fid("/home/john/Desktop/recon_native_debug.dat",std::ios::binary);
      fid.write((char*)m_reconstructed_data.get(),m_rp.nx * m_rp.ny * m_num_slices_native*sizeof(float));
    }

    // Free any data that is no longer needed
    cudaFree(m_d_filtered_projection_data);
    cudaFree(d_table_positions);
    cudaFree(d_tube_angles);
  }

  void ReconstructionMachine::ApplyFinalSliceThicknessing(){

    //// Copy data back to host
    //m_reconstructed_data.reset(new float[m_rp.nx * m_rp.ny * m_slice_locations_requested_slice_width.size()]);
    //cudaMemcpy(m_reconstructed_data.get(),d_reconstructed_slices_native_thickness)
  }

  void ReconstructionMachine::PrintCTGeometry(){
    std::cout << "CT Geometry and Scan derived parameters: "   << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "Num projections per turn:           "        << m_cg.projections_per_rotation     << std::endl;
    std::cout << "Num detector channels:              "        << m_cg.num_detector_cols            << std::endl;
    std::cout << "Num detector rows:                  "        << m_cg.num_detector_rows            << std::endl;
    std::cout << "Radius src->isocenter (mm):         "        << m_cg.distance_source_to_isocenter << std::endl;
    std::cout << "Radius src->detector (mm):          "        << m_cg.distance_source_to_detector  << std::endl;
    std::cout << "Table feed per rotation (mm):       "        << m_cg.z_rot                        << std::endl;
    std::cout << "Theta cone (rad):                   "        << m_cg.theta_cone                   << std::endl;
    std::cout << "Central channel:                    "        << m_cg.detector_central_col         << std::endl;
    std::cout << "Central row:                        "        << m_cg.detector_central_row         << std::endl;
    std::cout << "Acquisition FOV (mm):               "        << m_cg.acquisition_field_of_view    << std::endl;  
    std::cout << "Collimated slicewidth at isocenter: "        << m_cg.collimated_slice_width       << std::endl;
  }
}

namespace{
  cufftComplex * generate_filter(CTGeometry cg, float c, float a){
    // Create a spatial domain ramp filter.  Eventually we'll expost c and a
    // so users can customize filter response for smoother/sharper reconstructions  
    //float ds = mr->cg.r_f*sin(mr->cg.fan_angle_increment/2.0f); // This is at isocenter.  Is that correct?

    std::cout << "Generating filter" << std::endl;

    // Allocate the host filter
    std::shared_ptr<float> h_filter(new float[cg.num_detector_cols]);

    // Calculate the filter
    float pi_f = 3.141592653589f;
    float ds = cg.detector_pixel_size_col;
  
    auto r = [](float t)->float{
               float v = sin(t)/t + (cos(t)-1.0f)/(t*t);
               if (t==0)
                 v=0.5;
               return v;
             };

    int N = cg.num_detector_cols;
  
    for (int i=-N/2;i<N/2;i++){    
      h_filter.get()[i+N/2] = (c*c/(2.0f*ds)) * (a*r(c*pi_f*i) +
                                                 (((1.0f-a)/2.0f)*r(pi_f*c*i + pi_f)) +
                                                 (((1.0f -a)/2.0f)*r(pi_f*c*i-pi_f)));
    }

    // Apply the "fftshift" operation
    for (int i=0; i<N/2;i++){
      float tmp = h_filter.get()[i];   
      h_filter.get()[i]     = h_filter.get()[i+N/2];
      h_filter.get()[i+N/2] = tmp;
    }
  
    // Send to device
    cudaError_t cuda_status;
    cufftResult cufft_status;
  
    float * d_filter;
    cuda_status = cudaMalloc(&d_filter,cg.num_detector_cols*sizeof(float));
    gpuErrChk(cuda_status);

    cuda_status = cudaMemcpy(d_filter,h_filter.get(),cg.num_detector_cols*sizeof(float),cudaMemcpyHostToDevice);
    gpuErrChk(cuda_status);

    // Take the FFT to get it into the Fourier domain, and return pointer to the complex FFT array
    cufftHandle plan;
    cufft_status = cufftPlan1d(&plan,cg.num_detector_cols,CUFFT_R2C,1);
    cufftErrChk(cufft_status);
    
    cufftComplex * d_filter_final;
    int output_size = N/2+1;
  
    cuda_status = cudaMalloc(&d_filter_final,output_size*sizeof(cufftComplex));
    gpuErrChk(cuda_status);
  
    cufft_status = cufftExecR2C(plan,(cufftReal*)d_filter,d_filter_final);
    cufftErrChk(cufft_status);

    cudaFree(d_filter);
    cufftDestroy(plan);

    return d_filter_final;
  }
}