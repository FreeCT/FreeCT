/* CTBangBang is GPU and CPU CT reconstruction Software */
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

__constant__ struct recon_info d_ri;
__constant__ struct recon_params d_rp;

__global__ void reorder_block(float * output,float * input,int block_idx){

    int nx=d_rp.nx;
    int ny=d_rp.ny;

    int x_idx=threadIdx.x+blockIdx.x*blockDim.x;
    int y_idx=threadIdx.y+blockIdx.y*blockDim.y;
    int z_idx=threadIdx.z+blockIdx.z*blockDim.z;
    
    size_t block_offset=block_idx*nx*ny*d_ri.n_slices_block;

    size_t out_idx=z_idx*nx*ny+y_idx*nx+x_idx+block_offset;
    size_t in_idx=((d_ri.n_slices_block-1)-z_idx)*nx*ny+y_idx*nx+x_idx+block_offset;

    output[out_idx]=input[in_idx];
    
}

__global__ void thicken_slices(float * final_image_stack,float * reordered_image_stack,float * raw_recon_locations,int final_slice_idx,float slice_location){

    extern __shared__ float weights[];

    int nx=d_rp.nx;
    int ny=d_rp.ny;

    int n_raw_images=d_ri.n_slices_block*d_ri.n_blocks;

    int x_idx=threadIdx.x+blockIdx.x*blockDim.x;
    int y_idx=threadIdx.y+blockIdx.y*blockDim.y;

    // Each thread computes a weight term
    size_t w_idx=threadIdx.y*blockDim.x+threadIdx.x;
    
    if (w_idx<n_raw_images)
	weights[w_idx]=fmaxf(0.0f,1.0f-fabsf(raw_recon_locations[w_idx]-slice_location)/d_rp.slice_thickness);
	   
    __syncthreads();

    float sum_weights=0.0f;
    for (int i=0;i<n_raw_images;i++){
	sum_weights+=weights[i];
    }

    int out_idx=final_slice_idx*nx*ny+y_idx*nx+x_idx;
    for (int raw_slice=0;raw_slice<n_raw_images;raw_slice++){
	if (weights[raw_slice]!=0){
	    int raw_idx=raw_slice*nx*ny+y_idx*nx+x_idx;
	    final_image_stack[out_idx]+=(weights[raw_slice]/sum_weights)*reordered_image_stack[raw_idx];
	}
    }    
}
