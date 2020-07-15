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

#include <recon_structs.h>

#define pi 3.1415926535897f
#define N_PIX 2 // There may be a time in the future were this needs to be more flexible, but for right now, 2 works for all scanners

texture<float,cudaTextureType2D,cudaReadModeElementType> tex_a;
texture<float,cudaTextureType2D,cudaReadModeElementType> tex_b;
texture<float,cudaTextureType2D,cudaReadModeElementType> tex_c;
texture<float,cudaTextureType2D,cudaReadModeElementType> tex_d;

__constant__ struct ct_geom d_cg;
__constant__ struct recon_info d_ri;
__constant__ float d_filter[3680];

/* --- Helper functions (called by kernels) --- */
__device__ inline float angle(float x1,float x2,float y1,float y2){
    return asin((x1*y2-x2*y1)/(sqrt(x1*x1+x2*x2)*sqrt(y1*y1+y2*y2)));
}

__device__ inline float beta_rk(float da,float dr,float channel,int os_flag){
    float b0=(channel-pow(2.0f,os_flag)*d_cg.central_channel)*(d_cg.fan_angle_increment/pow(2.0f,os_flag));
    return angle(-(d_cg.r_f+dr),-(da),-(d_cg.src_to_det*cos(b0)+dr),-(d_cg.src_to_det*sin(b0)+da));
}

__device__ inline float d_alpha_r(float da,float dr){
    return angle(d_cg.r_f,0,d_cg.r_f+dr,da);
}

__device__ inline float r_fr(float da, float dr){
    return sqrt((d_cg.r_f+dr)*(d_cg.r_f+dr)+da*da);
}

__device__ inline float get_beta_idx(float beta,float * beta_lookup,int n_elements){
    int idx_low=0;

    while (beta>beta_lookup[idx_low]&&idx_low<(n_elements-1)){
    	idx_low++;
    }

    if (idx_low==0)
	idx_low++; 
    
    return (float)idx_low-1.0f+(beta-beta_lookup[idx_low-1])/(beta_lookup[idx_low]-beta_lookup[idx_low-1]);
}

__global__ void beta_lookup(float * lookup,float dr, float da,int os_flag){
    int channel=threadIdx.x+blockIdx.x*blockDim.x;   
    lookup[channel]=beta_rk(da,dr,channel,os_flag);
}

/* --- No flying focal spot rebinning kernels --- */
__global__ void n1_rebin(float * output,int row){
    int channel = threadIdx.x + blockDim.x*blockIdx.x;
    int proj    = threadIdx.y + blockDim.y*blockIdx.y;
    
    float beta=asin(((float)channel-2*d_cg.central_channel)*(d_cg.fan_angle_increment/2));    
    float alpha_idx=(float)proj-beta*d_cg.n_proj_turn/(2.0f*pi);
    float beta_idx=beta/d_cg.fan_angle_increment+d_cg.central_channel;

    int out_idx=d_cg.n_channels_oversampled*(d_ri.n_proj_pull/d_ri.n_ffs)*row+(d_ri.n_proj_pull/d_ri.n_ffs)*channel+proj;
    
    output[out_idx]=tex2D(tex_a,beta_idx+0.5f,alpha_idx+0.5f);
}

__global__ void n2_rebin(float * output,int row){
    int channel = threadIdx.x + blockDim.x*blockIdx.x;
    int proj    = threadIdx.y + blockDim.y*blockIdx.y;

    float beta=asin(((float)channel-2*d_cg.central_channel)*(d_cg.fan_angle_increment/2));
    float alpha_idx=(float)proj-beta*d_cg.n_proj_turn/(2.0f*pi);
    float beta_idx=beta/d_cg.fan_angle_increment+d_cg.central_channel;

    int out_idx=d_cg.n_channels_oversampled*(d_ri.n_proj_pull/d_ri.n_ffs)*row+(d_ri.n_proj_pull/d_ri.n_ffs)*channel+proj;
    
    output[out_idx]=tex2D(tex_b,beta_idx+0.5f,alpha_idx+0.5f);
}

/* --- Phi only flying focal spot rebinning kernels ---*/
__global__ void p_reshape(float * raw, float * out,int offset){
    int channel=blockIdx.x*blockDim.x+threadIdx.x;
    int row=blockIdx.y*blockDim.y+threadIdx.y;
    int proj=blockIdx.z;

    int idx_out=(proj+offset/2)+(d_ri.n_proj_pull/d_ri.n_ffs*channel)+(d_ri.n_proj_pull/d_ri.n_ffs*d_cg.n_channels*row);
    int idx_out_offset=idx_out+(d_cg.n_rows_raw*d_cg.n_channels*d_ri.n_proj_pull/d_ri.n_ffs);
    
    out[idx_out]       =raw[d_cg.n_channels*d_cg.n_rows_raw *(2*proj) +d_cg.n_channels*row+channel];
    out[idx_out_offset]=raw[d_cg.n_channels*d_cg.n_rows_raw*(2*proj+1)+d_cg.n_channels*row+channel];
}

__global__ void p1_rebin_t(float * output,float da,int row,float * beta_lookup){
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;
    int out_idx = d_cg.n_channels_oversampled*n_proj*row+n_proj*(2*channel)+proj;

    da=da;
   
    float beta = beta_rk(da,0,channel,0);
    beta_lookup[2*channel]=beta;
    float alpha_idx=d_ri.n_ffs*proj-beta*d_cg.n_proj_ffs/(2.0f*pi)-d_alpha_r(da,0)*d_cg.n_proj_ffs/(2.0f*pi);
    alpha_idx=alpha_idx/2.0f;
    
    output[out_idx]=tex2D(tex_a,alpha_idx+0.5f,channel+0.5f);
}

__global__ void p2_rebin_t(float * output,float da,int row,float * beta_lookup){
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;
    int out_idx = d_cg.n_channels_oversampled*n_proj*row+n_proj*(2*channel+1)+proj;

    da=-da;
    
    float beta = beta_rk(da,0,channel,0);
    beta_lookup[2*channel+1]=beta;
    float alpha_idx=d_ri.n_ffs*proj-beta*d_cg.n_proj_ffs/(2.0f*pi)-d_alpha_r(da,0)*d_cg.n_proj_ffs/(2.0f*pi);
    alpha_idx=(alpha_idx-1.0f)/2.0f;
    
    output[out_idx]=tex2D(tex_b,alpha_idx+0.5f,channel+0.5f);
}

__global__ void p1_rebin(float* output,float da,int row,float * beta_lookup){
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;
    int out_idx = d_cg.n_channels_oversampled*n_proj*row+n_proj*channel+proj;

    float beta  = asin((channel-2*d_cg.central_channel)*(d_cg.fan_angle_increment/2));
    float beta_idx=get_beta_idx(beta,beta_lookup,d_cg.n_channels_oversampled);
    
    output[out_idx]=tex2D(tex_a,proj+0.5f,beta_idx+0.5f); 
}

__global__ void p2_rebin(float* output,float da,int row,float * beta_lookup){
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;
    int out_idx = d_cg.n_channels_oversampled*n_proj*row+n_proj*channel+proj;

    float beta  = asin((channel-2*d_cg.central_channel)*(d_cg.fan_angle_increment/2));
    float beta_idx=get_beta_idx(beta,beta_lookup,d_cg.n_channels_oversampled);
    
    output[out_idx]=tex2D(tex_b,proj+0.5f,beta_idx+0.5f);     
}


/* --- Z only flying focal spot rebinning kernels ---*/
__global__ void z_reshape(float * raw, float * out,int offset){
    int channel=blockIdx.x*blockDim.x+threadIdx.x;
    int row=blockIdx.y*blockDim.y+threadIdx.y;
    int proj=blockIdx.z;

    int idx_out=(proj+offset/2)+(d_ri.n_proj_pull/d_ri.n_ffs*channel)+(d_ri.n_proj_pull/d_ri.n_ffs*d_cg.n_channels*row);
    int idx_out_offset=idx_out+(d_cg.n_rows_raw*d_cg.n_channels*d_ri.n_proj_pull/d_ri.n_ffs);
    
    out[idx_out]       =raw[d_cg.n_channels*d_cg.n_rows_raw *(2*proj) +d_cg.n_channels*row+channel];
    out[idx_out_offset]=raw[d_cg.n_channels*d_cg.n_rows_raw*(2*proj+1)+d_cg.n_channels*row+channel];
}


__global__ void z1_rebin(float * output,float * beta_lookup,float dr,int row){
    // This kernel handles all projections coming from focal spot 1
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    float beta=asin((channel-2.0f*d_cg.central_channel)*(d_cg.fan_angle_increment/2.0f)*d_cg.r_f/r_fr(0.0f,-dr));
    float alpha_idx=d_ri.n_ffs*proj-beta*d_cg.n_proj_ffs/(2.0f*pi)-d_alpha_r(0.0f,-dr)*d_cg.n_proj_ffs/(2.0f*pi);
    float beta_idx=get_beta_idx(beta,beta_lookup,d_cg.n_channels);

    alpha_idx=alpha_idx/2.0f;
    
    __syncthreads();

    int out_idx;
    if (!d_cg.reverse_row_interleave)
	out_idx=(d_ri.n_proj_pull/d_ri.n_ffs)*d_cg.n_channels_oversampled*2*row+(d_ri.n_proj_pull/d_ri.n_ffs)*channel+proj;
    else
	out_idx=(d_ri.n_proj_pull/d_ri.n_ffs)*d_cg.n_channels_oversampled*(2*row+1)+(d_ri.n_proj_pull/d_ri.n_ffs)*channel+proj;
    
    output[out_idx]=tex2D(tex_a,alpha_idx+0.5f,beta_idx+0.5f);
}

__global__ void z2_rebin(float * output,float * beta_lookup,float dr,int row){
    // This kernel handles all projections coming from focal spot 2
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    float beta=asin((channel-2.0f*d_cg.central_channel)*(d_cg.fan_angle_increment/2.0f)*d_cg.r_f/r_fr(0.0f,dr));
    float alpha_idx=d_ri.n_ffs*proj-beta*d_cg.n_proj_ffs/(2.0f*pi)-d_alpha_r(0.0f,dr)*d_cg.n_proj_ffs/(2.0f*pi);
    float beta_idx=get_beta_idx(beta,beta_lookup,d_cg.n_channels);

    alpha_idx=(alpha_idx-1.0f)/2.0f;
    
    __syncthreads();

    int out_idx;
    if (!d_cg.reverse_row_interleave)
	out_idx=(d_ri.n_proj_pull/d_ri.n_ffs)*d_cg.n_channels_oversampled*(2*row+1)+(d_ri.n_proj_pull/d_ri.n_ffs)*channel+proj;
    else
	out_idx=(d_ri.n_proj_pull/d_ri.n_ffs)*d_cg.n_channels_oversampled*2*row+(d_ri.n_proj_pull/d_ri.n_ffs)*channel+proj;    
    
    output[out_idx]=tex2D(tex_b,alpha_idx+0.5f,beta_idx+0.5f);
}

/* --- Z and Phi flying focal spot rebinning kernels ---*/
__global__ void a_reshape(float * raw, float * out,int offset){
    int channel=blockIdx.x*blockDim.x+threadIdx.x;
    int row=blockIdx.y*blockDim.y+threadIdx.y;
    int proj=blockIdx.z;

    int idx_out_1=(proj+offset/4)+(d_ri.n_proj_pull/d_ri.n_ffs*channel)+(d_ri.n_proj_pull/d_ri.n_ffs*d_cg.n_channels*row);
    int idx_out_2=idx_out_1+(d_cg.n_rows_raw*d_cg.n_channels*d_ri.n_proj_pull/d_ri.n_ffs);
    int idx_out_3=idx_out_2+(d_cg.n_rows_raw*d_cg.n_channels*d_ri.n_proj_pull/d_ri.n_ffs);
    int idx_out_4=idx_out_3+(d_cg.n_rows_raw*d_cg.n_channels*d_ri.n_proj_pull/d_ri.n_ffs);    
    
    out[idx_out_1]=raw[d_cg.n_channels*d_cg.n_rows_raw *(4*proj) +d_cg.n_channels*row+channel];
    out[idx_out_2]=raw[d_cg.n_channels*d_cg.n_rows_raw*(4*proj+1)+d_cg.n_channels*row+channel];
    out[idx_out_3]=raw[d_cg.n_channels*d_cg.n_rows_raw*(4*proj+2)+d_cg.n_channels*row+channel];
    out[idx_out_4]=raw[d_cg.n_channels*d_cg.n_rows_raw*(4*proj+3)+d_cg.n_channels*row+channel];
}

__global__ void a1_rebin_t(float * output,float da, float dr, int row,float * beta_lookup){
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;
    int out_idx = d_cg.n_channels_oversampled*n_proj*row+n_proj*(2*channel)+proj;

    float beta=beta_rk(da,-dr,channel,0);
    beta_lookup[2*channel]=beta;
    //float alpha_idx=(proj)-beta*d_cg.n_proj_turn/(2.0f*pi)-d_alpha_r(da,-dr)*d_cg.n_proj_turn/(2.0f*pi);
    float alpha_idx=d_ri.n_ffs*proj-beta*d_cg.n_proj_ffs/(2.0f*pi)-d_alpha_r(da,-dr)*d_cg.n_proj_ffs/(2.0f*pi);
    alpha_idx=alpha_idx/4.0f;
    
    output[out_idx]=tex2D(tex_a,alpha_idx+0.5f,channel+0.5f);
}

__global__ void a2_rebin_t(float * output,float da, float dr, int row,float * beta_lookup){
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;
    int out_idx = d_cg.n_channels_oversampled*n_proj*row+n_proj*(2*channel+1)+proj;

    float beta=beta_rk(-da,-dr,channel,0);
    beta_lookup[2*channel+1]=beta;
    //float alpha_idx=(proj)-beta*d_cg.n_proj_turn/(2.0f*pi)-d_alpha_r(-da,-dr)*d_cg.n_proj_turn/(2.0f*pi);
    float alpha_idx=d_ri.n_ffs*proj-beta*d_cg.n_proj_ffs/(2.0f*pi)-d_alpha_r(-da,-dr)*d_cg.n_proj_ffs/(2.0f*pi);
    alpha_idx=(alpha_idx-1.0f)/4.0f;
    
    output[out_idx]=tex2D(tex_b,alpha_idx+0.5f,channel+0.5f);
}

__global__ void a3_rebin_t(float * output,float da, float dr, int row,float * beta_lookup){
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;
    int out_idx = d_cg.n_channels_oversampled*n_proj*row+n_proj*(2*channel)+proj;

    float beta=beta_rk(da,dr,channel,0);
    beta_lookup[2*channel]=beta;    
    //float alpha_idx=(proj)-beta*d_cg.n_proj_turn/(2.0f*pi)-d_alpha_r(da,dr)*d_cg.n_proj_turn/(2.0f*pi);
    float alpha_idx=d_ri.n_ffs*proj-beta*d_cg.n_proj_ffs/(2.0f*pi)-d_alpha_r(da,dr)*d_cg.n_proj_ffs/(2.0f*pi);
    alpha_idx=(alpha_idx-2.0f)/4.0f;

    output[out_idx]=tex2D(tex_c,alpha_idx+0.5f,channel+0.5f);
}

__global__ void a4_rebin_t(float * output,float da, float dr, int row,float * beta_lookup){
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;
    int out_idx = d_cg.n_channels_oversampled*n_proj*row+n_proj*(2*channel+1)+proj;

    float beta=beta_rk(-da,dr,channel,0);
    beta_lookup[2*channel+1]=beta;
    //float alpha_idx=(proj)-beta*d_cg.n_proj_turn/(2.0f*pi)-d_alpha_r(-da,dr)*d_cg.n_proj_turn/(2.0f*pi);
    float alpha_idx=d_ri.n_ffs*proj-beta*d_cg.n_proj_ffs/(2.0f*pi)-d_alpha_r(-da,dr)*d_cg.n_proj_ffs/(2.0f*pi);
    alpha_idx=(alpha_idx-3.0f)/4.0f;
    
    output[out_idx]=tex2D(tex_d,alpha_idx+0.5f,channel+0.5f);
}

__global__ void a1_rebin_b(float * output,float * beta_lookup,float dr,int row){
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;

    int out_idx;
    if (!d_cg.reverse_row_interleave)
	out_idx = n_proj*d_cg.n_channels_oversampled*2*row+n_proj*channel+proj;
    else
	out_idx = n_proj*d_cg.n_channels_oversampled*(2*row+1)+n_proj*channel+proj;
   
    float beta=asin((channel-2.0f*d_cg.central_channel)*(d_cg.fan_angle_increment/2.0f)*d_cg.r_f/r_fr(0.0f,-dr));

    float beta_idx=get_beta_idx(beta,beta_lookup,d_cg.n_channels_oversampled);
    
    __syncthreads();

    output[out_idx]=tex2D(tex_a,proj+0.5f,beta_idx+0.5f);
}

__global__ void a2_rebin_b(float * output,float * beta_lookup,float dr,int row){
    int channel = threadIdx.x+blockIdx.x*blockDim.x;
    int proj    = threadIdx.y+blockIdx.y*blockDim.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;

    int out_idx;
    if (!d_cg.reverse_row_interleave)
	out_idx = n_proj*d_cg.n_channels_oversampled*(2*row+1)+n_proj*channel+proj;
    else
	out_idx = n_proj*d_cg.n_channels_oversampled*2*row+n_proj*channel+proj;
    
    float beta=asin((channel-2.0f*d_cg.central_channel)*(d_cg.fan_angle_increment/2.0f)*d_cg.r_f/r_fr(0.0f,dr));
    float beta_idx=get_beta_idx(beta,beta_lookup,d_cg.n_channels_oversampled);
    
    __syncthreads();

    output[out_idx]=tex2D(tex_b,proj+0.5f,beta_idx+0.5f);
}

/* --- Reshape out kernel ---*/
// Reshapes row-sheet rebinned array into projection-sheet array
__global__ void reshape_out(float * output,float * input){

    size_t offset=d_cg.add_projections;

    size_t k = threadIdx.x+blockIdx.x*blockDim.x;//j channel
    size_t j = threadIdx.y+blockIdx.y*blockDim.y;//k proj 
    size_t i = threadIdx.z+blockIdx.z*blockDim.z;//i row

    size_t in_idx=(d_cg.n_channels_oversampled*d_ri.n_proj_pull/d_ri.n_ffs)*i+d_ri.n_proj_pull/d_ri.n_ffs*j+(k+offset);
    size_t out_idx=k*d_cg.n_channels_oversampled*d_cg.n_rows+i*d_cg.n_channels_oversampled+j;

    output[out_idx]=input[in_idx];
}

/* --- Filter kernel --- */
__global__ void filter(float * output, int row){

    extern __shared__ float s[];
    float * s_row=s;
    float * s_filter=(float*) &s_row[d_cg.n_channels_oversampled];

    float conv_pix[N_PIX]={0};
    
    int channel=threadIdx.x*N_PIX;
    int proj=blockIdx.y;

    int n_proj  = d_ri.n_proj_pull/d_ri.n_ffs;
    int first_output_pixel=row*d_cg.n_channels_oversampled*n_proj+n_proj*channel+proj;

    int N=d_cg.n_channels_oversampled;
    
    // Load in 2*n_pix of filter data and n_pix of row data to shared memory
#pragma unroll
    for (int i=0;i<N_PIX;i++){
	s_row[channel+i]=output[first_output_pixel+i*n_proj];
	s_filter[channel+i]=d_filter[channel+i];
	s_filter[channel+i+N]=d_filter[channel+i+N];
    }

    __syncthreads(); // Make sure every thread has finished copying data into shared memory

    // Compute n_pix of the convolution    
#pragma unroll
    for (int i=0;i<N_PIX;i++){
	int l=channel+i;
	for (int k=0;k<N;k++){
	    conv_pix[i]+=s_filter[l-k+(N)]*s_row[k];//s_filter[l-k+(N-1)]*
	}	
    }

    // Copy back to global memory
#pragma unroll
    for (int i=0;i<N_PIX;i++){
	output[first_output_pixel+i*n_proj]=conv_pix[i];
    }
}

