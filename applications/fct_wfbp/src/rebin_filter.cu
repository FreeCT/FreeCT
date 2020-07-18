/* FreeCT_wFBP is GPU and CPU CT reconstruction Software */
/* Copyright (C) 2015  John Hoffman */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

/* Questions and comments should be directed to */
/* jmhoffman@mednet.ucla.edu with "CTBANGBANG" in the subject line*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <ctbb_macros.h>
#include <recon_structs.h>
#include <rebin_filter.cuh>
#include <rebin_filter.h>

#include <iostream>
#include <fstream>

void copy_sheet(float * sheetptr, int row, struct recon_metadata *mr, struct ct_geom cg);
void load_filter(float * f_array,struct recon_metadata * mr);
void generate_filter(float * f_array,struct recon_metadata * mr, float c=1.0f, float a=1.0f); // c and a control ramp rolloff. Not exposed at this point.

void rebin_nffs(struct recon_metadata *mr);
void rebin_pffs(struct recon_metadata *mr);
void rebin_zffs(struct recon_metadata *mr);
void rebin_affs(struct recon_metadata *mr);

int rebin_filter(struct recon_metadata * mr){

  std::cout << "REBIN FILTER: " << mr->ri.n_ffs << std::endl;
  
    switch (mr->ri.n_ffs){
    case 1:{
	rebin_nffs(mr);
	break;}
    case 2:{
	if (mr->rp.z_ffs==1)
	    rebin_zffs(mr);
	else
	    rebin_pffs(mr);
	break;}
    case 4:{
	rebin_affs(mr);
	break;}
    }
    
    return 0;
}

void rebin_nffs(struct recon_metadata *mr){

    if (mr->flags.timing){
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start);
    }

    cudaStream_t stream1,stream2;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);
    
    struct ct_geom cg=mr->cg;
    struct recon_info ri=mr->ri;
    
    // Allocate the entire output array here and on GPU Note that size
    // of output array is the same for all FFS configurations
    float * d_output;

    cudaMalloc(&d_output,cg.n_channels_oversampled*cg.n_rows*mr->ri.n_proj_pull/mr->ri.n_ffs*sizeof(float));
    cudaMemset(d_output,0,cg.n_channels_oversampled*cg.n_rows*mr->ri.n_proj_pull/mr->ri.n_ffs*sizeof(float));
    
    // Copy of ct geometry
    cudaMemcpyToSymbol(d_cg,&cg,sizeof(struct ct_geom),0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_ri,&ri,sizeof(struct recon_info),0,cudaMemcpyHostToDevice);
    
    // Ready our filter
    float * h_filter=(float*)calloc(2*cg.n_channels_oversampled,sizeof(float));
    //load_filter(h_filter,mr);
    generate_filter(h_filter,mr);

    //std::ofstream filter_fid("filter_test.bin",std::ios::binary);
    //filter_fid.write((char*)h_filter,2*cg.n_channels_oversampled*sizeof(float));

    cudaMemcpyToSymbol(d_filter,h_filter,2*cg.n_channels_oversampled*sizeof(float),0,cudaMemcpyHostToDevice);

    // Configure textures (see rebin_filter.cuh)
    tex_a.addressMode[0] = cudaAddressModeClamp;
    tex_a.addressMode[1] = cudaAddressModeClamp;
    tex_a.addressMode[2] = cudaAddressModeClamp;
    tex_a.filterMode     = cudaFilterModeLinear;
    tex_a.normalized     = false;

    tex_b.addressMode[0] = cudaAddressModeClamp;
    tex_b.addressMode[1] = cudaAddressModeClamp;
    tex_b.addressMode[2] = cudaAddressModeClamp;
    tex_b.filterMode     = cudaFilterModeLinear;
    tex_b.normalized     = false;

    // Projection data will be thought of as row-sheets (since we row-wise rebin)
    size_t proj_array_size=cg.n_channels*mr->ri.n_proj_pull;
    cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc<float>();

    // Allocate raw data array 1
    float * sheet_1=(float *)calloc(cg.n_channels*mr->ri.n_proj_pull,sizeof(float));
    cudaArray * cu_raw_1;
    cudaMallocArray(&cu_raw_1,&channelDesc,cg.n_channels,mr->ri.n_proj_pull);
    // Allocate raw data array 2
    float * sheet_2=(float *)calloc(cg.n_channels*mr->ri.n_proj_pull,sizeof(float));
    cudaArray * cu_raw_2;
    cudaMallocArray(&cu_raw_2,&channelDesc,cg.n_channels,mr->ri.n_proj_pull);

    // Kernel Dimensions
    dim3 rebin_threads(8,8);
    dim3 rebin_blocks(cg.n_channels_oversampled/rebin_threads.x,mr->ri.n_proj_pull/rebin_threads.y);

    // Determine the number of threads needed based on the #define N_PIX (found in rebin_filter.cuh)
    int n_threads=cg.n_channels_oversampled/N_PIX;
    
    dim3 filter_threads(n_threads,1,1);
    dim3 filter_blocks(1,mr->ri.n_proj_pull/mr->ri.n_ffs,1);
    unsigned int shared_size=cg.n_channels_oversampled*sizeof(float)+2*cg.n_channels_oversampled*sizeof(float); // row+filter;
    
    // Reshape raw data into row sheets
    std::cout << cg.n_channels << std::endl;
    std::cout << cg.n_rows_raw << std::endl;
    std::cout << cg.n_rows << std::endl;
    std::cout << mr->ri.n_proj_pull << std::endl;
    float * sheets=(float*)calloc(cg.n_channels*cg.n_rows_raw*mr->ri.n_proj_pull,sizeof(float));
    for (int i=0;i<(int)cg.n_rows_raw;i++){
      for (int j=0;j<(int)cg.n_channels;j++){
        for (int k=0;k<(int)mr->ri.n_proj_pull;k++){
                int in_idx = k*(int)cg.n_channels*(int)cg.n_rows_raw+i*(int)cg.n_channels+j;
                int out_idx= i*(int)cg.n_channels*(int)mr->ri.n_proj_pull+k*(int)cg.n_channels+((int)cg.n_channels - 1 -j);
                sheets[out_idx]=mr->ctd.raw[in_idx];
	    }
	}
    }
	
    for (int i=0;i<cg.n_rows;i+=2){
	// Copy first set of projections over to GPU	
	cudaMemcpyToArrayAsync(cu_raw_1,0,0,&sheets[i*proj_array_size],proj_array_size*sizeof(float),cudaMemcpyHostToDevice,stream1);
	cudaBindTextureToArray(tex_a,cu_raw_1,channelDesc);

	// Launch Kernel A
	n1_rebin<<<rebin_blocks,rebin_threads,0,stream1>>>(d_output,i);
	filter<<<filter_blocks,filter_threads,shared_size,stream1>>>(d_output,i);
	    
	//Begin the second transfer while 1st kernel executes
	cudaMemcpyToArrayAsync(cu_raw_2,0,0,&sheets[(i+1)*proj_array_size],proj_array_size*sizeof(float),cudaMemcpyHostToDevice,stream2);
	cudaBindTextureToArray(tex_b,cu_raw_2,channelDesc);

	n2_rebin<<<rebin_blocks,rebin_threads,0,stream2>>>(d_output,i+1);
	filter<<<filter_blocks,filter_threads,shared_size,stream2>>>(d_output,i+1);
    }
	
    cudaFreeArray(cu_raw_1);
    cudaFreeArray(cu_raw_2);

    //Reshape data into our mr structure
    cudaMalloc(&mr->ctd.d_rebin,cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections)*sizeof(float));

    n_threads=1;
    size_t n_proj_out=(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections);
    dim3 threads_reshape_out(n_threads,1,1);
    dim3 blocks_reshape_out(n_proj_out/n_threads,cg.n_channels_oversampled,cg.n_rows);

    reshape_out<<<blocks_reshape_out,threads_reshape_out>>>(mr->ctd.d_rebin,d_output);

    // Check "testing" flag, write rebin to disk if set
    //if (mr->flags.testing){
    //cudaMemcpy(mr->ctd.rebin,mr->ctd.d_rebin,cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections)*sizeof(float),cudaMemcpyDeviceToHost);
    //char fullpath[4096+255];
    //strcpy(fullpath,mr->rp.output_dir);
    //strcat(fullpath,"/rebin.ct_test");
    //FILE * outfile=fopen(fullpath,"w");
    //fwrite(mr->ctd.rebin,sizeof(float),cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull-2*cg.add_projections_ffs)/mr->ri.n_ffs,outfile);
    //fclose(outfile);
    //}

    cudaFree(d_output);
    
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);

    free(h_filter);
}

void rebin_pffs(struct recon_metadata *mr){

    // Set up some constants on the host 
    struct ct_geom cg=mr->cg;
    struct recon_info ri=mr->ri;

    const double da=cg.src_to_det*cg.r_f*cg.fan_angle_increment/(4.0f*(cg.src_to_det-cg.r_f));
    
    // Set up some constants and infrastructure on the GPU
    cudaStream_t stream1,stream2;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);

    cudaMemcpyToSymbol(d_cg,&cg,sizeof(struct ct_geom),0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_ri,&ri,sizeof(struct recon_info),0,cudaMemcpyHostToDevice);
    
    int proj_per_call=32;

    // Need to split arrays by focal spot and reshape into "sheets"
    float * d_raw_1;
    float * d_raw_2;
    gpuErrchk(cudaMalloc(&d_raw_1,proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float)));
    gpuErrchk(cudaMalloc(&d_raw_2,proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float)));

    float * d_fs;
    gpuErrchk(cudaMalloc(&d_fs,mr->ri.n_proj_pull*cg.n_channels*cg.n_rows_raw*sizeof(float)));
    cudaMemset(d_fs,0,mr->ri.n_proj_pull*cg.n_channels*cg.n_rows_raw*sizeof(float));
    
    //dim3 threads_reshape(32,16,1);
    dim3 threads_reshape(8,16,1);
    dim3 blocks_reshape(cg.n_channels/threads_reshape.x,cg.n_rows_raw/threads_reshape.y,proj_per_call/2);
    
    int i=0;
    while (i<mr->ri.n_proj_pull){

	cudaMemcpyAsync(d_raw_1,&mr->ctd.raw[i*cg.n_channels*cg.n_rows_raw],proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float),cudaMemcpyHostToDevice,stream1);
	p_reshape<<<blocks_reshape,threads_reshape,0,stream1>>>(d_raw_1,d_fs,i);
	
	cudaMemcpyAsync(d_raw_2,&mr->ctd.raw[(i+proj_per_call)*cg.n_channels*cg.n_rows_raw],proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float),cudaMemcpyHostToDevice,stream2);
	p_reshape<<<blocks_reshape,threads_reshape,0,stream2>>>(d_raw_2,d_fs,i+proj_per_call);

	i+=2*proj_per_call;
    }

    if (mr->flags.testing){
	float * h_fs=(float*)calloc(mr->ri.n_proj_pull*cg.n_rows_raw*cg.n_channels,sizeof(float));
	cudaMemcpy(h_fs,d_fs,mr->ri.n_proj_pull*cg.n_rows_raw*cg.n_channels*sizeof(float),cudaMemcpyDeviceToHost);
	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/h_fs.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(h_fs,sizeof(float),mr->ri.n_proj_pull*cg.n_rows_raw*cg.n_channels,outfile);
	fclose(outfile);
	free(h_fs);
    }
    
    cudaFree(d_raw_1);
    cudaFree(d_raw_2);

    // Configure textures (see rebin_filter.cuh)
    tex_a.addressMode[0] = cudaAddressModeClamp;
    tex_a.addressMode[1] = cudaAddressModeClamp;
    tex_a.addressMode[2] = cudaAddressModeClamp;
    tex_a.filterMode     = cudaFilterModeLinear;
    tex_a.normalized     = false;

    tex_b.addressMode[0] = cudaAddressModeClamp;
    tex_b.addressMode[1] = cudaAddressModeClamp;
    tex_b.addressMode[2] = cudaAddressModeClamp;
    tex_b.filterMode     = cudaFilterModeLinear;
    tex_b.normalized     = false;

    cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc<float>();

    cudaArray * cu_raw_1;
    cudaArray * cu_raw_2;
    gpuErrchk(cudaMallocArray(&cu_raw_1,&channelDesc,ri.n_proj_pull/ri.n_ffs,cg.n_channels));
    gpuErrchk(cudaMallocArray(&cu_raw_2,&channelDesc,ri.n_proj_pull/ri.n_ffs,cg.n_channels));
    
    float * d_rebin_t;
    gpuErrchk(cudaMalloc(&d_rebin_t,cg.n_channels_oversampled*cg.n_rows*ri.n_proj_pull/ri.n_ffs*sizeof(float)));
    cudaMemset(d_rebin_t,0,cg.n_channels_oversampled*cg.n_rows*ri.n_proj_pull/ri.n_ffs*sizeof(float));
    
    //dim3 threads_t_rebin(32,32);
    dim3 threads_t_rebin(8,32);    
    dim3 blocks_t_rebin(cg.n_channels/threads_t_rebin.x,ri.n_proj_pull/ri.n_ffs/threads_t_rebin.y);

    // Set up our lookup table for use in the channel rebinning
    float * d_beta_lookup;
    gpuErrchk(cudaMalloc(&d_beta_lookup,cg.n_channels_oversampled*sizeof(float)));

    for (int i=0;i<cg.n_rows;i++){
	cudaMemcpyToArrayAsync(cu_raw_1,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream1);
	cudaBindTextureToArray(tex_a,cu_raw_1,channelDesc);

	cudaMemcpyToArrayAsync(cu_raw_2,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i+cg.n_channels*ri.n_proj_pull/ri.n_ffs*cg.n_rows_raw],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream2);
	cudaBindTextureToArray(tex_b,cu_raw_2,channelDesc);

	p1_rebin_t<<<blocks_t_rebin,threads_t_rebin,0,stream1>>>(d_rebin_t,da,i,d_beta_lookup);
	p2_rebin_t<<<blocks_t_rebin,threads_t_rebin,0,stream2>>>(d_rebin_t,da,i,d_beta_lookup);
     }
    
    cudaFree(d_fs);
    cudaFreeArray(cu_raw_1);
    cudaFreeArray(cu_raw_2);
    
    if (mr->flags.testing){
	float * h_rebin_t=(float*)calloc(cg.n_channels_oversampled*cg.n_rows*ri.n_proj_pull/ri.n_ffs,sizeof(float));
	cudaMemcpy(h_rebin_t,d_rebin_t,cg.n_channels_oversampled*cg.n_rows*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToHost);
	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/rebin_t.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(h_rebin_t,sizeof(float),cg.n_channels_oversampled*cg.n_rows*ri.n_proj_pull/ri.n_ffs,outfile);
	fclose(outfile);
	free(h_rebin_t);
    }

    gpuErrchk(cudaMallocArray(&cu_raw_1,&channelDesc,ri.n_proj_pull/ri.n_ffs,cg.n_channels_oversampled));
    gpuErrchk(cudaMallocArray(&cu_raw_2,&channelDesc,ri.n_proj_pull/ri.n_ffs,cg.n_channels_oversampled));

    float * d_output;
    gpuErrchk(cudaMalloc(&d_output,cg.n_channels_oversampled*cg.n_rows*ri.n_proj_pull/ri.n_ffs*sizeof(float)));

    // Ready our filter
    float * h_filter=(float*)calloc(2*cg.n_channels_oversampled,sizeof(float));

    load_filter(h_filter,mr);
    gpuErrchk(cudaMemcpyToSymbol(d_filter,h_filter,2*cg.n_channels_oversampled*sizeof(float),0,cudaMemcpyHostToDevice));
	
    int sheet_size=cg.n_channels_oversampled*ri.n_proj_pull/ri.n_ffs;

    //dim3 threads_rebin(32,32);
    dim3 threads_rebin(16,32);
    dim3 blocks_rebin(cg.n_channels_oversampled/threads_rebin.x,ri.n_proj_pull/ri.n_ffs/threads_rebin.y);
    
    // Determine the number of threads needed based on the #define N_PIX (found in rebin_filter.cuh)
    int n_threads=cg.n_channels_oversampled/N_PIX;
    
    dim3 threads_filter(n_threads,1,1);
    dim3 blocks_filter(1,mr->ri.n_proj_pull/mr->ri.n_ffs,1);
    unsigned int shared_size=cg.n_channels_oversampled*sizeof(float)+2*cg.n_channels_oversampled*sizeof(float); // row+filter;

    for (int i=0;i<cg.n_rows;i+=2){

	gpuErrchk(cudaMemcpyToArrayAsync(cu_raw_1,0,0,&d_rebin_t[i*sheet_size],sheet_size*sizeof(float),cudaMemcpyDeviceToDevice,stream1));
	gpuErrchk(cudaBindTextureToArray(tex_a,cu_raw_1,channelDesc));

	gpuErrchk(cudaMemcpyToArrayAsync(cu_raw_2,0,0,&d_rebin_t[(i+1)*sheet_size],sheet_size*sizeof(float),cudaMemcpyDeviceToDevice,stream2));
	gpuErrchk(cudaBindTextureToArray(tex_b,cu_raw_2,channelDesc));

	p1_rebin<<<blocks_rebin,threads_rebin,0,stream1>>>(d_output,da,i,d_beta_lookup);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	filter<<<blocks_filter,threads_filter,shared_size,stream1>>>(d_output,i);
	    
	p2_rebin<<<blocks_rebin,threads_rebin,0,stream2>>>(d_output,da,i+1,d_beta_lookup);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	filter<<<blocks_filter,threads_filter,shared_size,stream2>>>(d_output,i+1);
	
    }

    cudaFree(d_rebin_t);

    cudaMalloc(&mr->ctd.d_rebin,cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections)*sizeof(float));

    n_threads=128;
    size_t n_proj_out=(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections);
    dim3 threads_reshape_out(n_threads,1,1);
    dim3 blocks_reshape_out(n_proj_out/n_threads,cg.n_channels_oversampled,cg.n_rows);        

    reshape_out<<<blocks_reshape_out,threads_reshape_out>>>(mr->ctd.d_rebin,d_output);
    


    // Check "testing" flag, write rebin to disk if set
    if (mr->flags.testing){
	cudaMemcpy(mr->ctd.rebin,mr->ctd.d_rebin,cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections)*sizeof(float),cudaMemcpyDeviceToHost);
	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/rebin.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(mr->ctd.rebin,sizeof(float),cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull-2*cg.add_projections_ffs)/mr->ri.n_ffs,outfile);
	fclose(outfile);
    }

    cudaFree(d_beta_lookup);
    cudaFree(d_output);
    cudaFreeArray(cu_raw_1);
    cudaFreeArray(cu_raw_2);
    
    free(h_filter);
    
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    
}
void rebin_zffs(struct recon_metadata *mr){


    // Two GPU parts to this rebin:
    // (1) Split arrays and reshape into sheets
    // (2) Rebin sheets and interleave into final array

    struct ct_geom cg=mr->cg;
    struct recon_info ri=mr->ri;
    struct recon_params rp=mr->rp;

    const double da=0.0;
    const double dr=cg.src_to_det*rp.coll_slicewidth/(4.0*(cg.src_to_det-cg.r_f)*tan(cg.anode_angle));
    
    cudaStream_t stream1,stream2;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);

    // Copy of ct geometry
    cudaMemcpyToSymbol(d_cg,&cg,sizeof(struct ct_geom),0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_ri,&ri,sizeof(struct recon_info),0,cudaMemcpyHostToDevice);
    
    int proj_per_call=32;

    // Need to split arrays by focal spot and reshape into "sheets"
    float * d_raw_1;
    float * d_raw_2;
    cudaMalloc(&d_raw_1,proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float));
    cudaMalloc(&d_raw_2,proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float));
    
    float * d_fs;
    cudaMalloc(&d_fs,mr->ri.n_proj_pull*cg.n_channels*cg.n_rows_raw*sizeof(float));

    dim3 threads_reshape(32,16,1);
    dim3 blocks_reshape(cg.n_channels/threads_reshape.x,cg.n_rows_raw/threads_reshape.y,proj_per_call/2);
    
    int i=0;
    while (i<mr->ri.n_proj_pull){
	cudaMemcpyAsync(d_raw_1,&mr->ctd.raw[i*cg.n_channels*cg.n_rows_raw],proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float),cudaMemcpyHostToDevice,stream1);
	z_reshape<<<blocks_reshape,threads_reshape,0,stream1>>>(d_raw_1,d_fs,i);
	
	cudaMemcpyAsync(d_raw_2,&mr->ctd.raw[(i+proj_per_call)*cg.n_channels*cg.n_rows_raw],proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float),cudaMemcpyHostToDevice,stream2);
	z_reshape<<<blocks_reshape,threads_reshape,0,stream2>>>(d_raw_2,d_fs,i+proj_per_call);

	i+=2*proj_per_call;
    }

    cudaFree(d_raw_1);
    cudaFree(d_raw_2);

    // Check "testing" flag, write rebin to disk if set
    if (mr->flags.testing){
	float * h_fs_1;
	float * h_fs_2;
	h_fs_1=(float*)calloc((mr->ri.n_proj_pull/mr->ri.n_ffs)*cg.n_channels*cg.n_rows_raw,sizeof(float));
	h_fs_2=(float*)calloc((mr->ri.n_proj_pull/mr->ri.n_ffs)*cg.n_channels*cg.n_rows_raw,sizeof(float));

	cudaMemcpyAsync(h_fs_1, d_fs                                                     ,cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToHost,stream1);
	cudaMemcpyAsync(h_fs_2,&d_fs[cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs],cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToHost,stream2);    

	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/reshape_1.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(h_fs_1,sizeof(float),cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs,outfile);
	fclose(outfile);

	memset(fullpath,0,4096+255);
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/reshape_2.ct_test");
	outfile=fopen(fullpath,"w");
	fwrite(h_fs_1,sizeof(float),cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs,outfile);
	fclose(outfile);

	free(h_fs_1);
	free(h_fs_2);
    }

    
    // Step 1 finished
    // Set up everything for the rebinning
    // Configure textures (see rebin_filter.cuh)
    tex_a.addressMode[0] = cudaAddressModeClamp;
    tex_a.addressMode[1] = cudaAddressModeClamp;
    tex_a.addressMode[2] = cudaAddressModeClamp;
    tex_a.filterMode     = cudaFilterModeLinear;
    tex_a.normalized     = false;

    tex_b.addressMode[0] = cudaAddressModeClamp;
    tex_b.addressMode[1] = cudaAddressModeClamp;
    tex_b.addressMode[2] = cudaAddressModeClamp;
    tex_b.filterMode     = cudaFilterModeLinear;
    tex_b.normalized     = false;

    // Projection data will be thought of as row-sheets (since we row-wise rebin)
    //size_t proj_array_size=cg.n_channels*mr->ri.n_proj_pull;
    cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc<float>();

    cudaArray * cu_raw_1;
    cudaArray * cu_raw_2;
    cudaMallocArray(&cu_raw_1,&channelDesc,mr->ri.n_proj_pull/ri.n_ffs,cg.n_channels);
    cudaMallocArray(&cu_raw_2,&channelDesc,mr->ri.n_proj_pull/ri.n_ffs,cg.n_channels);

    float * d_output;
    cudaMalloc(&d_output,cg.n_channels_oversampled*cg.n_rows*ri.n_proj_pull/ri.n_ffs*sizeof(float));

    // Allocate and compute beta lookup tables
    float * d_beta_lookup_1;
    float * d_beta_lookup_2;
    cudaMalloc(&d_beta_lookup_1,cg.n_channels*sizeof(float));
    cudaMalloc(&d_beta_lookup_2,cg.n_channels*sizeof(float));
    beta_lookup<<<1,cg.n_channels>>>(d_beta_lookup_1,-dr,da,0);
    beta_lookup<<<1,cg.n_channels>>>(d_beta_lookup_2, dr,da,0);

    if (mr->flags.testing){
	float * h_b_lookup_1=(float*)calloc(cg.n_channels,sizeof(float));
	float * h_b_lookup_2=(float*)calloc(cg.n_channels,sizeof(float));
	
	cudaMemcpy(h_b_lookup_1,d_beta_lookup_1,cg.n_channels*sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemcpy(h_b_lookup_2,d_beta_lookup_2,cg.n_channels*sizeof(float),cudaMemcpyDeviceToHost);

	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/beta_lookup_1.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(h_b_lookup_1,sizeof(float),cg.n_channels,outfile);
	fwrite(h_b_lookup_2,sizeof(float),cg.n_channels,outfile);
	fclose(outfile);

	free(h_b_lookup_1);
	free(h_b_lookup_2);
    }

    // Ready our filter
    float * h_filter=(float*)calloc(2*cg.n_channels_oversampled,sizeof(float));
    load_filter(h_filter,mr);
    cudaMemcpyToSymbol(d_filter,h_filter,2*cg.n_channels_oversampled*sizeof(float),0,cudaMemcpyHostToDevice);
    
    dim3 threads_rebin(32,32);
    dim3 blocks_rebin(cg.n_channels_oversampled/threads_rebin.x,ri.n_proj_pull/ri.n_ffs/threads_rebin.y);

    // Determine the number of threads needed based on the #define N_PIX (found in rebin_filter.cuh)
    int n_threads=cg.n_channels_oversampled/N_PIX;
    
    dim3 threads_filter(n_threads,1,1);
    dim3 blocks_filter(1,mr->ri.n_proj_pull/mr->ri.n_ffs,1);
    unsigned int shared_size=cg.n_channels_oversampled*sizeof(float)+2*cg.n_channels_oversampled*sizeof(float); // row+filter;
    
    for (int i=0;i<cg.n_rows_raw;i++){
	cudaMemcpyToArrayAsync(cu_raw_1,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream1);
	cudaBindTextureToArray(tex_a,cu_raw_1,channelDesc);

	cudaMemcpyToArrayAsync(cu_raw_2,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i+cg.n_channels*ri.n_proj_pull/ri.n_ffs*cg.n_rows_raw],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream2);
	cudaBindTextureToArray(tex_b,cu_raw_2,channelDesc);

	z1_rebin<<<blocks_rebin,threads_rebin,0,stream1>>>(d_output,d_beta_lookup_1,dr,i);
	z2_rebin<<<blocks_rebin,threads_rebin,0,stream2>>>(d_output,d_beta_lookup_2,dr,i);

	cudaDeviceSynchronize();

	filter<<<blocks_filter,threads_filter,shared_size,stream1>>>(d_output,2*i);
	filter<<<blocks_filter,threads_filter,shared_size,stream2>>>(d_output,2*i+1);
    }

    //Reshape into array ready for backprojection
    cudaFree(d_fs);
    cudaMalloc(&mr->ctd.d_rebin,cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections)*sizeof(float));

    n_threads=128;
    size_t n_proj_out=(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections);
    dim3 threads_reshape_out(n_threads,1,1);
    dim3 blocks_reshape_out(n_proj_out/n_threads,cg.n_channels_oversampled,cg.n_rows);        

    reshape_out<<<blocks_reshape_out,threads_reshape_out>>>(mr->ctd.d_rebin,d_output);

    // Check "testing" flag, write rebin to disk if set
    if (mr->flags.testing){
	cudaMemcpy(mr->ctd.rebin,mr->ctd.d_rebin,cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections)*sizeof(float),cudaMemcpyDeviceToHost);
	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/rebin.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(mr->ctd.rebin,sizeof(float),cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull-2*cg.add_projections_ffs)/mr->ri.n_ffs,outfile);
	fclose(outfile);
    }

    //free(h_output);
    free(h_filter);
    //cudaFree(d_fs);
    cudaFree(d_output);
    cudaFree(d_beta_lookup_1);
    cudaFree(d_beta_lookup_2);
    
    cudaFreeArray(cu_raw_1);
    cudaFreeArray(cu_raw_2);
    
    cudaStreamDestroy(stream1); 
    cudaStreamDestroy(stream2);
    
}
void rebin_affs(struct recon_metadata *mr){

    // Set up some constants on the host 
    struct ct_geom cg=mr->cg;
    struct recon_info ri=mr->ri;
    struct recon_params rp=mr->rp;

    const double da=cg.src_to_det*cg.r_f*cg.fan_angle_increment/(4.0f*(cg.src_to_det-cg.r_f));
    const double dr=cg.src_to_det*rp.coll_slicewidth/(4.0*(cg.src_to_det-cg.r_f)*tan(cg.anode_angle));
    
    // Set up some constants and infrastructure on the GPU
    cudaStream_t stream1,stream2,stream3,stream4;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);
    cudaStreamCreate(&stream3);
    cudaStreamCreate(&stream4);

    cudaMemcpyToSymbol(d_cg,&cg,sizeof(struct ct_geom),0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_ri,&ri,sizeof(struct recon_info),0,cudaMemcpyHostToDevice);

    // Configure textures (see rebin_filter.cuh)
    tex_a.addressMode[0] = cudaAddressModeClamp;
    tex_a.addressMode[1] = cudaAddressModeClamp;
    tex_a.addressMode[2] = cudaAddressModeClamp;
    tex_a.filterMode     = cudaFilterModeLinear;
    tex_a.normalized     = false;

    tex_b.addressMode[0] = cudaAddressModeClamp;
    tex_b.addressMode[1] = cudaAddressModeClamp;
    tex_b.addressMode[2] = cudaAddressModeClamp;
    tex_b.filterMode     = cudaFilterModeLinear;
    tex_b.normalized     = false;

    tex_c.addressMode[0] = cudaAddressModeClamp;
    tex_c.addressMode[1] = cudaAddressModeClamp;
    tex_c.addressMode[2] = cudaAddressModeClamp;
    tex_c.filterMode     = cudaFilterModeLinear;
    tex_c.normalized     = false;

    tex_d.addressMode[0] = cudaAddressModeClamp;
    tex_d.addressMode[1] = cudaAddressModeClamp;
    tex_d.addressMode[2] = cudaAddressModeClamp;
    tex_d.filterMode     = cudaFilterModeLinear;
    tex_d.normalized     = false;

    cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc<float>();

    // Need to split arrays by focal spot and reshape into "sheets"
    int proj_per_call=32;

    float * d_raw_1;
    float * d_raw_2;
    float * d_raw_3;
    float * d_raw_4;
    cudaMalloc(&d_raw_1,proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float));
    cudaMalloc(&d_raw_2,proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float));
    cudaMalloc(&d_raw_3,proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float));
    cudaMalloc(&d_raw_4,proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float));
    
    float * d_fs;
    cudaMalloc(&d_fs,mr->ri.n_proj_pull*cg.n_channels*cg.n_rows_raw*sizeof(float));

    dim3 threads_reshape(32,16,1);
    dim3 blocks_reshape(cg.n_channels/threads_reshape.x,cg.n_rows_raw/threads_reshape.y,proj_per_call/4);
    
    int i=0;
    while (i<mr->ri.n_proj_pull){

	cudaMemcpyAsync(d_raw_1,&mr->ctd.raw[i*cg.n_channels*cg.n_rows_raw],proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float),cudaMemcpyHostToDevice,stream1);
	a_reshape<<<blocks_reshape,threads_reshape,0,stream1>>>(d_raw_1,d_fs,i);
	
	cudaMemcpyAsync(d_raw_2,&mr->ctd.raw[(i+proj_per_call)*cg.n_channels*cg.n_rows_raw],proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float),cudaMemcpyHostToDevice,stream2);
	a_reshape<<<blocks_reshape,threads_reshape,0,stream2>>>(d_raw_2,d_fs,i+proj_per_call);

	cudaMemcpyAsync(d_raw_3,&mr->ctd.raw[(i+2*proj_per_call)*cg.n_channels*cg.n_rows_raw],proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float),cudaMemcpyHostToDevice,stream3);
	a_reshape<<<blocks_reshape,threads_reshape,0,stream3>>>(d_raw_3,d_fs,i+2*proj_per_call);

	cudaMemcpyAsync(d_raw_4,&mr->ctd.raw[(i+3*proj_per_call)*cg.n_channels*cg.n_rows_raw],proj_per_call*cg.n_channels*cg.n_rows_raw*sizeof(float),cudaMemcpyHostToDevice,stream4);
	a_reshape<<<blocks_reshape,threads_reshape,0,stream4>>>(d_raw_4,d_fs,i+3*proj_per_call);

	i+=4*proj_per_call;
    }

    // Check "testing" flag, write rebin to disk if set
    if (mr->flags.testing){
	float * h_fs_1;
	float * h_fs_2;
	float * h_fs_3;
	float * h_fs_4;
	h_fs_1=(float*)calloc((mr->ri.n_proj_pull/mr->ri.n_ffs)*cg.n_channels*cg.n_rows_raw,sizeof(float));
	h_fs_2=(float*)calloc((mr->ri.n_proj_pull/mr->ri.n_ffs)*cg.n_channels*cg.n_rows_raw,sizeof(float));
	h_fs_3=(float*)calloc((mr->ri.n_proj_pull/mr->ri.n_ffs)*cg.n_channels*cg.n_rows_raw,sizeof(float));
	h_fs_4=(float*)calloc((mr->ri.n_proj_pull/mr->ri.n_ffs)*cg.n_channels*cg.n_rows_raw,sizeof(float));

	cudaMemcpyAsync(h_fs_1, d_fs                                                     ,cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToHost,stream1);
	cudaMemcpyAsync(h_fs_2,&d_fs[  cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs],cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToHost,stream2);    
	cudaMemcpyAsync(h_fs_3,&d_fs[2*cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs],cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToHost,stream3);
	cudaMemcpyAsync(h_fs_4,&d_fs[3*cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs],cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToHost,stream4);    

	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/reshape_1.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(h_fs_1,sizeof(float),cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs,outfile);
	fclose(outfile);

	memset(fullpath,0,4096+255);
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/reshape_2.ct_test");
	outfile=fopen(fullpath,"w");
	fwrite(h_fs_2,sizeof(float),cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs,outfile);
	fclose(outfile);

	memset(fullpath,0,4096+255);
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/reshape_3.ct_test");
	outfile=fopen(fullpath,"w");
	fwrite(h_fs_3,sizeof(float),cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs,outfile);
	fclose(outfile);

	memset(fullpath,0,4096+255);
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/reshape_4.ct_test");
	outfile=fopen(fullpath,"w");
	fwrite(h_fs_4,sizeof(float),cg.n_channels*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs,outfile);
	fclose(outfile);

	free(h_fs_1);
	free(h_fs_2);
	free(h_fs_3);
	free(h_fs_4);
    }
    
    cudaFree(d_raw_1);
    cudaFree(d_raw_2);
    cudaFree(d_raw_3);
    cudaFree(d_raw_4);

    float * d_rebin_t1;
    float * d_rebin_t2;
    cudaMalloc(&d_rebin_t1,ri.n_proj_pull/ri.n_ffs*cg.n_channels_oversampled*cg.n_rows_raw*sizeof(float));
    cudaMalloc(&d_rebin_t2,ri.n_proj_pull/ri.n_ffs*cg.n_channels_oversampled*cg.n_rows_raw*sizeof(float));

    cudaArray * cu_raw_1;
    cudaArray * cu_raw_2;
    cudaArray * cu_raw_3;
    cudaArray * cu_raw_4;    
    cudaMallocArray(&cu_raw_1,&channelDesc,ri.n_proj_pull/ri.n_ffs,cg.n_channels);
    cudaMallocArray(&cu_raw_2,&channelDesc,ri.n_proj_pull/ri.n_ffs,cg.n_channels);
    cudaMallocArray(&cu_raw_3,&channelDesc,ri.n_proj_pull/ri.n_ffs,cg.n_channels);
    cudaMallocArray(&cu_raw_4,&channelDesc,ri.n_proj_pull/ri.n_ffs,cg.n_channels);

    dim3 threads_t_rebin(32,32);
    dim3 blocks_t_rebin(cg.n_channels/threads_t_rebin.x,ri.n_proj_pull/ri.n_ffs/threads_t_rebin.y);

    // Allocate and compute beta lookup tables for channel rebinning
    float * d_beta_lookup_1;
    float * d_beta_lookup_2;
    cudaMalloc(&d_beta_lookup_1,cg.n_channels_oversampled*sizeof(float));
    cudaMalloc(&d_beta_lookup_2,cg.n_channels_oversampled*sizeof(float));
    
    for (int i=0;i<cg.n_rows_raw;i++){
	cudaMemcpyToArrayAsync(cu_raw_1,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i                                                        ],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream1);
	cudaMemcpyToArrayAsync(cu_raw_2,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i +   cg.n_channels*ri.n_proj_pull/ri.n_ffs*cg.n_rows_raw],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream2);
	cudaMemcpyToArrayAsync(cu_raw_3,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i + 2*cg.n_channels*ri.n_proj_pull/ri.n_ffs*cg.n_rows_raw],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream3);
	cudaMemcpyToArrayAsync(cu_raw_4,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i + 3*cg.n_channels*ri.n_proj_pull/ri.n_ffs*cg.n_rows_raw],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream4);

	cudaBindTextureToArray(tex_a,cu_raw_1,channelDesc);
	a1_rebin_t<<<blocks_t_rebin,threads_t_rebin,0,stream1>>>(d_rebin_t1,da,dr,i,d_beta_lookup_1);
	
	cudaMemcpyToArrayAsync(cu_raw_2,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i+cg.n_channels*ri.n_proj_pull/ri.n_ffs*cg.n_rows_raw],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream2);
	cudaBindTextureToArray(tex_b,cu_raw_2,channelDesc);
	a2_rebin_t<<<blocks_t_rebin,threads_t_rebin,0,stream2>>>(d_rebin_t1,da,dr,i,d_beta_lookup_1);
	
	cudaMemcpyToArrayAsync(cu_raw_3,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i+2*cg.n_channels*ri.n_proj_pull/ri.n_ffs*cg.n_rows_raw],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream3);
	cudaBindTextureToArray(tex_c,cu_raw_3,channelDesc);
	a3_rebin_t<<<blocks_t_rebin,threads_t_rebin,0,stream3>>>(d_rebin_t2,da,dr,i,d_beta_lookup_2);
	
	cudaMemcpyToArrayAsync(cu_raw_4,0,0,&d_fs[cg.n_channels*ri.n_proj_pull/ri.n_ffs*i+3*cg.n_channels*ri.n_proj_pull/ri.n_ffs*cg.n_rows_raw],cg.n_channels*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream4);
	cudaBindTextureToArray(tex_d,cu_raw_4,channelDesc);
	a4_rebin_t<<<blocks_t_rebin,threads_t_rebin,0,stream4>>>(d_rebin_t2,da,dr,i,d_beta_lookup_2);
	
    }

    if (mr->flags.testing){
	float * h_rebin_t1;
	float * h_rebin_t2;
	h_rebin_t1=(float*)calloc(cg.n_channels_oversampled*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs,sizeof(float));
	h_rebin_t2=(float*)calloc(cg.n_channels_oversampled*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs,sizeof(float));
	
	cudaMemcpy(h_rebin_t1,d_rebin_t1,cg.n_channels_oversampled*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemcpy(h_rebin_t2,d_rebin_t2,cg.n_channels_oversampled*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToHost);

	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/rebin_t1.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(h_rebin_t1,sizeof(float),cg.n_channels_oversampled*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs,outfile);
	fclose(outfile);
	
	memset(fullpath,0,4096+255);
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/rebin_t2.ct_test");
	outfile=fopen(fullpath,"w");
	fwrite(h_rebin_t2,sizeof(float),cg.n_channels_oversampled*cg.n_rows_raw*ri.n_proj_pull/ri.n_ffs,outfile);
	fclose(outfile);

	free(h_rebin_t1);
	free(h_rebin_t2);	
    }

    
    cudaFree(d_fs);
    cudaFreeArray(cu_raw_1);
    cudaFreeArray(cu_raw_2);
    cudaFreeArray(cu_raw_3);
    cudaFreeArray(cu_raw_4);

    float * d_output;
    gpuErrchk(cudaMalloc(&d_output,cg.n_channels_oversampled*cg.n_rows*ri.n_proj_pull/ri.n_ffs*sizeof(float)));

    gpuErrchk(cudaMallocArray(&cu_raw_1,&channelDesc,ri.n_proj_pull/ri.n_ffs,cg.n_channels_oversampled));
    gpuErrchk(cudaMallocArray(&cu_raw_2,&channelDesc,ri.n_proj_pull/ri.n_ffs,cg.n_channels_oversampled));

    if (mr->flags.testing){
	float * h_b_lookup_1=(float*)calloc(cg.n_channels_oversampled,sizeof(float));
	float * h_b_lookup_2=(float*)calloc(cg.n_channels_oversampled,sizeof(float));
	
	cudaMemcpy(h_b_lookup_1,d_beta_lookup_1,cg.n_channels_oversampled*sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemcpy(h_b_lookup_2,d_beta_lookup_2,cg.n_channels_oversampled*sizeof(float),cudaMemcpyDeviceToHost);
	
	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/beta_lookup_1.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(h_b_lookup_1,sizeof(float),cg.n_channels_oversampled,outfile);
	fclose(outfile);

	memset(fullpath,0,4096+255);
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/beta_lookup_2.ct_test");
	fopen(fullpath,"w");
	fwrite(h_b_lookup_2,sizeof(float),cg.n_channels_oversampled,outfile);
	fclose(outfile);
	
	free(h_b_lookup_1);
	free(h_b_lookup_2);
    }

    // Ready our filter
    float * h_filter=(float*)calloc(2*cg.n_channels_oversampled,sizeof(float));
    load_filter(h_filter,mr);
    cudaMemcpyToSymbol(d_filter,h_filter,2*cg.n_channels_oversampled*sizeof(float),0,cudaMemcpyHostToDevice);
    
    if (mr->flags.testing){
	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/filter.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(h_filter,sizeof(float),2*cg.n_channels_oversampled,outfile);
	fclose(outfile);
    }
    
    dim3 threads_rebin(32,32);
    dim3 blocks_rebin(cg.n_channels_oversampled/threads_rebin.x,ri.n_proj_pull/ri.n_ffs/threads_rebin.y);

    // Determine the number of threads needed based on the #define N_PIX (found in rebin_filter.cuh)
    int n_threads=cg.n_channels_oversampled/N_PIX;
    
    dim3 threads_filter(n_threads,1,1);
    dim3 blocks_filter(1,mr->ri.n_proj_pull/mr->ri.n_ffs,1);
    unsigned int shared_size=cg.n_channels_oversampled*sizeof(float)+2*cg.n_channels_oversampled*sizeof(float); // row+filter;
    
    for (int i=0;i<cg.n_rows_raw;i++){
	cudaMemcpyToArrayAsync(cu_raw_1,0,0,&d_rebin_t1[cg.n_channels_oversampled*ri.n_proj_pull/ri.n_ffs*i],cg.n_channels_oversampled*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream1);
	cudaBindTextureToArray(tex_a,cu_raw_1,channelDesc);

	cudaMemcpyToArrayAsync(cu_raw_2,0,0,&d_rebin_t2[cg.n_channels_oversampled*ri.n_proj_pull/ri.n_ffs*i],cg.n_channels_oversampled*ri.n_proj_pull/ri.n_ffs*sizeof(float),cudaMemcpyDeviceToDevice,stream2);
	cudaBindTextureToArray(tex_b,cu_raw_2,channelDesc);

	a1_rebin_b<<<blocks_rebin,threads_rebin,0,stream1>>>(d_output,d_beta_lookup_1,dr,i);	
	a2_rebin_b<<<blocks_rebin,threads_rebin,0,stream2>>>(d_output,d_beta_lookup_2,dr,i);

	cudaDeviceSynchronize();
	
	filter<<<blocks_filter,threads_filter,shared_size,stream1>>>(d_output,2*i);
	filter<<<blocks_filter,threads_filter,shared_size,stream2>>>(d_output,2*i+1);
    }

    //Reshape data into our mr structure
    cudaFree(d_rebin_t1);
    cudaFree(d_rebin_t2);
    
    cudaMalloc(&mr->ctd.d_rebin,cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections)*sizeof(float));
    
    n_threads=128;
    size_t n_proj_out=(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections);
    dim3 threads_reshape_out(n_threads,1,1);
    dim3 blocks_reshape_out(n_proj_out/n_threads,cg.n_channels_oversampled,cg.n_rows);
    
    reshape_out<<<blocks_reshape_out,threads_reshape_out>>>(mr->ctd.d_rebin,d_output);
    
    // Check "testing" flag, write rebin to disk if set
    if (mr->flags.testing){
	cudaMemcpy(mr->ctd.rebin,mr->ctd.d_rebin,cg.n_channels_oversampled*cg.n_rows*(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections)*sizeof(float),cudaMemcpyDeviceToHost);
	char fullpath[4096+255];
	strcpy(fullpath,mr->rp.output_dir);
	strcat(fullpath,"/rebin.ct_test");
	FILE * outfile=fopen(fullpath,"w");
	fwrite(mr->ctd.rebin,sizeof(float),cg.n_channels_oversampled*cg.n_rows*(ri.n_proj_pull-2*cg.add_projections_ffs)/ri.n_ffs,outfile);
	fclose(outfile);
    }

    free(h_filter);

    cudaFreeArray(cu_raw_1);
    cudaFreeArray(cu_raw_2);
    
    cudaFree(d_output);
    cudaFree(d_beta_lookup_1);
    cudaFree(d_beta_lookup_2);

    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    cudaStreamDestroy(stream4);
}

void copy_sheet(float * sheetptr, int row,struct recon_metadata * mr,struct ct_geom cg){
    for (int j=0;j<cg.n_channels;j++){
	for (int k=0;k<mr->ri.n_proj_pull;k++){
	    sheetptr[j+k*cg.n_channels]=mr->ctd.raw[k*cg.n_channels*cg.n_rows_raw+row*cg.n_channels+j];
	}
    }
}

void load_filter(float * f_array,struct recon_metadata * mr){
    struct ct_geom cg=mr->cg;
    struct recon_params rp=mr->rp;

    char fullpath[4096+255]={0};

    FILE * filter_file;
    switch (rp.recon_kernel){
    case -100:{
	sprintf(fullpath,"%s/resources/filters/f_%lu_ramp.txt",mr->install_dir,cg.n_channels);
	break;}
    case -1:{
	sprintf(fullpath,"%s/resources/filters/f_%lu_exp.txt",mr->install_dir,cg.n_channels);
	break;}
    case 1:{
	sprintf(fullpath,"%s/resources/filters/f_%lu_smooth.txt",mr->install_dir,cg.n_channels);
	break;}
    case 2:{
	sprintf(fullpath,"%s/resources/filters/f_%lu_medium.txt",mr->install_dir,cg.n_channels);
	break;}
    case 3:{
	sprintf(fullpath,"%s/resources/filters/f_%lu_sharp.txt",mr->install_dir,cg.n_channels);
	break;}
    default:{
	sprintf(fullpath,"%s/resources/filters/f_%lu_b%i.txt",mr->install_dir,cg.n_channels,rp.recon_kernel);
	break;}
    }

    filter_file=fopen(fullpath,"r");

    if (filter_file==NULL){
	perror("Filter file not found.  You may need to create one.  We hope to add on-the-fly creation in a future update.");
	exit(1);
    }

    fread(f_array,sizeof(float),2*cg.n_channels_oversampled,filter_file);
    fclose(filter_file);
}


void generate_filter(float * f_array,struct recon_metadata * mr, float c, float a){
  // Create a spatial domain ramp filter.  Eventually we'll expost c and a
  // so users can customize filter response for smoother/sharper reconstructions
  
  //float * h_filter=(float*)calloc(2*cg.n_channels_oversampled,sizeof(float));
  //float ds = mr->cg.r_f*sin(mr->cg.fan_angle_increment/2.0f); // This is at isocenter.  Is that correct?

  float pi_f = 3.141592653589f;
  float ds = mr->cg.src_to_det*sin(mr->cg.fan_angle_increment/2.0f); // I think it should be at the detector

  auto r = [](float t)-> float{
             float v = sin(t)/t + (cos(t)-1.0f)/(t*t);
             if (t==0)
               v=0.5;
             return v;
           };
  
  int test = (int)mr->cg.n_channels_oversampled;
  
  for (int i = -test;i < test;i++){
    

    f_array[i+mr->cg.n_channels_oversampled] = (c*c/(2.0f*ds)) * (a*r(c*pi_f*i) +
                                                                  (((1.0f-a)/2.0f)*r(pi_f*c*i + pi_f)) +
                                                                  (((1.0f -a)/2.0f)*r(pi_f*c*i-pi_f)));
    
  }
  
}
