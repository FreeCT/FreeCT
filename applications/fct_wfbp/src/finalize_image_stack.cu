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

#include <stdio.h>
#include <ctbb_macros.h>
#include <finalize_image_stack.cuh>
#include <finalize_image_stack.h>

int finalize_image_stack(struct recon_metadata * mr){

    struct recon_params rp=mr->rp;
    struct recon_info ri=mr->ri;

    cudaStream_t stream;
    cudaStreamCreate(&stream);
    
    // Copy reference structures to device
    cudaMemcpyToSymbol(d_ri,&mr->ri,sizeof(struct recon_info),0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_rp,&mr->rp,sizeof(struct recon_params),0,cudaMemcpyHostToDevice);

    float * d_raw_image_stack;
    cudaMalloc(&d_raw_image_stack,mr->rp.nx*mr->rp.ny*mr->ri.n_slices_recon*sizeof(float));
    cudaMemcpy(d_raw_image_stack,mr->ctd.image,mr->rp.nx*mr->rp.ny*mr->ri.n_slices_recon*sizeof(float),cudaMemcpyHostToDevice);

    float * d_temp_out;
    cudaMalloc(&d_temp_out,rp.nx*rp.ny*ri.n_slices_recon*sizeof(float));

    int recon_direction=fabs(rp.end_pos-rp.start_pos)/(rp.end_pos-rp.start_pos);
    if (recon_direction!=1&&recon_direction!=-1) // user request one slice (end_pos==start_pos)
	recon_direction=1;
    
    int table_direction=(mr->table_positions[1000]-mr->table_positions[0])/abs(mr->table_positions[1000]-mr->table_positions[0]);
    if (table_direction!=1&&table_direction!=-1)
	printf("Axial scans are currently unsupported, or a different error has occurred\n");
    
    // Check for a reversed stack of images and flip, otherwise just copy
    dim3 threads_reshape(1,1,ri.n_slices_block);
    dim3 blocks_reshape(rp.nx/threads_reshape.x,rp.ny/threads_reshape.y,1);
    if (recon_direction!=table_direction){
	for (int b=0;b<ri.n_blocks;b++){
	    reorder_block<<<blocks_reshape,threads_reshape>>>(d_temp_out,d_raw_image_stack,b);
	}
    }
    else{
	cudaMemcpy(d_temp_out,d_raw_image_stack,rp.nx*rp.ny*ri.n_slices_recon*sizeof(float),cudaMemcpyDeviceToDevice);
    }

    //float * temp_out=(float*)calloc(rp.nx*rp.ny*ri.n_slices_recon,sizeof(float));
    //cudaMemcpy(temp_out,d_temp_out,rp.nx*rp.ny*ri.n_slices_recon*sizeof(float),cudaMemcpyDeviceToHost);
    //FLOAT_DEBUG(temp_out,rp.nx*rp.ny*ri.n_slices_recon,"/home/john/Desktop/raw_array.bin");

    cudaFree(d_raw_image_stack);
    
    // Once we have straightened our image stack out, we need to adjust slice thickness
    // to match what the user requested.
    // We use a triangle average with the FWHM equal to the requested slice thickness

    int n_raw_images=ri.n_slices_block*ri.n_blocks;
    int n_slices_final=floor(fabs(rp.end_pos-rp.start_pos)/rp.slice_thickness)+1;
    mr->ctd.final_image_stack=(float*)calloc(rp.nx*rp.ny*n_slices_final,sizeof(float));
    cudaMalloc(&mr->ctd.d_final_image_stack,rp.nx*rp.ny*n_slices_final*sizeof(float));
    cudaMemset(mr->ctd.d_final_image_stack,0,rp.nx*rp.ny*n_slices_final*sizeof(float));
    
    float * recon_locations;
    recon_locations=(float*)calloc(n_slices_final,sizeof(float));
    for (int i=0;i<n_slices_final;i++){
	recon_locations[i]=rp.start_pos+recon_direction*i*rp.slice_thickness;
    }

    float * raw_recon_locations;
    float * d_raw_recon_locations;
    raw_recon_locations=(float*)calloc(n_raw_images,sizeof(float));
    cudaMalloc(&d_raw_recon_locations,n_raw_images*sizeof(float));
    for (int i=0;i<n_raw_images;i++){
	raw_recon_locations[i]=ri.recon_start_pos+recon_direction*i*rp.coll_slicewidth;//(rp.start_pos-recon_direction*rp.slice_thickness)+recon_direction*i*rp.coll_slicewidth;
    }
    cudaMemcpy(d_raw_recon_locations,raw_recon_locations,n_raw_images*sizeof(float),cudaMemcpyHostToDevice);

    //dim3 threads_thicken(32,32,1);
    dim3 threads_thicken(32,32,1);
    dim3 blocks_thicken(rp.nx/threads_thicken.x,rp.ny/threads_thicken.y,1);
    
    for (int k=0; k<n_slices_final; k++){
	float slice_location=recon_locations[k];
	thicken_slices<<<blocks_thicken,threads_thicken,n_raw_images*sizeof(float),stream>>>(mr->ctd.d_final_image_stack,d_temp_out,d_raw_recon_locations,k,slice_location);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
    }

    cudaMemcpy(mr->ctd.final_image_stack,mr->ctd.d_final_image_stack,rp.nx*rp.ny*n_slices_final*sizeof(float),cudaMemcpyDeviceToHost);

    // Free stuff allocated inside cleanup
    cudaFree(d_raw_recon_locations);
    cudaFree(d_temp_out);
    cudaFree(mr->ctd.d_final_image_stack);

    cudaStreamDestroy(stream);

    free(recon_locations);
    free(raw_recon_locations);
    
    return 0;
}
