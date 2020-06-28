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
#include <recon_structs.h>
#include <backproject.cuh>
#include <backproject.h>
#include <math.h>

#define Bx 16
#define By 16

int backproject(struct recon_metadata * mr){

    struct ct_geom cg=mr->cg;

    //float tube_start=mr->tube_angles[mr->ri.idx_pull_start+cg.add_projections_ffs]*pi/180;

    float tube_start=fmod((double)(mr->tube_angles[0]+((mr->ri.idx_pull_start+cg.add_projections_ffs)*360.0f/cg.n_proj_ffs)),360.0)*pi/180.0f;

    int n_half_turns=(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections)/(cg.n_proj_turn/2);
    
    cudaStream_t stream1,stream2;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);

    // Allocate the final volume array
    float * d_output;
    cudaMalloc(&d_output,mr->rp.nx*mr->rp.ny*mr->ri.n_slices_recon*sizeof(float));
    cudaMemset(d_output,0,mr->rp.nx*mr->rp.ny*mr->ri.n_slices_recon*sizeof(float));
    
    // Copy reference structures to device
    cudaMemcpyToSymbol(d_cg,&mr->cg,sizeof(struct ct_geom),0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_rp,&mr->rp,sizeof(struct recon_params),0,cudaMemcpyHostToDevice);
    
    // Configure textures (see backproject.cuh)
    tex_a.addressMode[0] = cudaAddressModeBorder;
    tex_a.addressMode[1] = cudaAddressModeBorder;
    tex_a.addressMode[2] = cudaAddressModeBorder;
    tex_a.filterMode     = cudaFilterModeLinear;
    tex_a.normalized     = false;

    tex_b.addressMode[0] = cudaAddressModeBorder;
    tex_b.addressMode[1] = cudaAddressModeBorder;
    tex_b.addressMode[2] = cudaAddressModeBorder;
    tex_b.filterMode     = cudaFilterModeLinear;
    tex_b.normalized     = false;

    cudaChannelFormatDesc channelDesc=cudaCreateChannelDesc<float>();

    cudaArray * cu_proj_1;
    cudaArray * cu_proj_2;
    cudaMallocArray(&cu_proj_1,&channelDesc,cg.n_channels_oversampled,I*cg.n_rows*n_half_turns);
    cudaMallocArray(&cu_proj_2,&channelDesc,cg.n_channels_oversampled,I*cg.n_rows*n_half_turns);

    dim3 threads(Bx,By);
    dim3 blocks(mr->rp.nx/Bx,mr->rp.ny/By,mr->ri.n_slices_block/K);

    for (int i=0;i<cg.n_proj_turn/2;i+=I*2){
	for (int k=0;k<n_half_turns;k++){
	    cudaMemcpyToArrayAsync(cu_proj_1,0,k*I*cg.n_rows,&mr->ctd.d_rebin[(i+k*cg.n_proj_turn/2)*cg.n_rows*cg.n_channels_oversampled],I*cg.n_rows*cg.n_channels_oversampled*sizeof(float),cudaMemcpyDeviceToDevice,stream1);
	}
	cudaBindTextureToArray(tex_a,cu_proj_1,channelDesc);

	for (int k=0;k<n_half_turns;k++){
	    cudaMemcpyToArrayAsync(cu_proj_2,0,k*I*cg.n_rows,&mr->ctd.d_rebin[(i+I+k*cg.n_proj_turn/2)*cg.n_rows*cg.n_channels_oversampled],I*cg.n_rows*cg.n_channels_oversampled*sizeof(float),cudaMemcpyDeviceToDevice,stream2);
	}
	cudaBindTextureToArray(tex_b,cu_proj_2,channelDesc);
	
	// Kernel call 1
	bp_a<<<blocks,threads,0,stream1>>>(d_output,i,tube_start,n_half_turns);

	// Kernel call 2
	bp_b<<<blocks,threads,0,stream2>>>(d_output,i+I,tube_start,n_half_turns);
	
    }

    long block_offset=(mr->ri.cb.block_idx-1)*mr->rp.nx*mr->rp.ny*mr->ri.n_slices_block;
    cudaMemcpy(&mr->ctd.image[block_offset],d_output,mr->rp.nx*mr->rp.ny*mr->ri.n_slices_block*sizeof(float),cudaMemcpyDeviceToHost);

    cudaFree(mr->ctd.d_rebin);
    cudaFree(d_output);
    cudaFreeArray(cu_proj_1);
    cudaFreeArray(cu_proj_2);
  
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);

    return 0;
}
