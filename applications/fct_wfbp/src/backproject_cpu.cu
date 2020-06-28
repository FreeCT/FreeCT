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

#include <interp.h>
#include <recon_structs.h>
#include <backproject_cpu.h>

#define pi 3.1415368979f

#define K 1
#define I 1 // Must be set to 1 currently

float W(float q){
    float out;
    float Q=0.6f;
    if (fabs(q)<Q){
	out=1.0f;
    }
    else if ((fabs(q)>=Q)&&(fabs(q)<1.0f)){
	out=pow(cos((pi/2.0f)*(fabs(q)-Q)/(1.0f-Q)),2.0f);
    }
    else {
	out=0.0f;
    }
    return out;
}

int backproject_cpu(struct recon_metadata * mr){
    
    struct ct_geom cg=mr->cg;
    struct recon_params rp=mr->rp;
    
    //float tube_start=mr->tube_angles[mr->ri.idx_pull_start+cg.add_projections_ffs]*pi/180;
    float tube_start=fmod((double)(mr->tube_angles[0]+((mr->ri.idx_pull_start+cg.add_projections_ffs)*360.0f/cg.n_proj_ffs)),360.0)*pi/180.0f;
    int n_half_turns=(mr->ri.n_proj_pull/mr->ri.n_ffs-2*cg.add_projections)/(cg.n_proj_turn/2);
    
    // Allocate the final output volume and intermediate voxel/weight arrays
    float * h_output;
    h_output=(float *)calloc(mr->rp.nx*mr->rp.ny*mr->ri.n_slices_block,sizeof(float));

    float * s_t;
    s_t=(float *)malloc(mr->rp.nx*mr->rp.ny*mr->ri.n_slices_block*sizeof(float));
    float * h_t;
    h_t=(float *)malloc(mr->rp.nx*mr->rp.ny*mr->ri.n_slices_block*sizeof(float));
    
    // Allocate the array to hold the projections currently being processed
    float * proj;
    proj=(float *)malloc(cg.n_channels_oversampled*cg.n_rows*I*n_half_turns*sizeof(float));

    // Info needed for interpolation from proj array
    struct array_dims dim;
    dim.idx1=cg.n_channels_oversampled;
    dim.idx2=cg.n_rows*n_half_turns;

    // For WFBP, outer loop is over angles in a half turn
    for (int i=0;i<cg.n_proj_turn/2;i+=I){

	printf("Backprojection projection set %d/%d...\n",i,cg.n_proj_turn/2);
	
	//Fetch n_half_turns projections for current projection angle
	for (int kk=0;kk<n_half_turns;kk++){
	    memcpy(&proj[cg.n_channels_oversampled*cg.n_rows*kk],&mr->ctd.rebin[(i+kk*cg.n_proj_turn/2)*cg.n_channels_oversampled*cg.n_rows],cg.n_channels_oversampled*cg.n_rows*sizeof(float));
	}
	
	// Set our intermediate arrays to zeros
	memset(s_t,0,mr->rp.nx*mr->rp.ny*mr->ri.n_slices_block*sizeof(float));
	memset(h_t,0,mr->rp.nx*mr->rp.ny*mr->ri.n_slices_block*sizeof(float));

	// Start loop over the voxels.
	// X and Y are looped over first to reuse theta, phat and p_idx calculations
	// as well as *theoretically* improve cache hit rates for the interpolations
	for (int xi=0;xi<mr->rp.nx;xi++){
	    float x=(rp.recon_fov/rp.nx)*((float)xi-(rp.nx-1)/2.0f)+rp.x_origin;

	    for (int yi=0;yi<mr->rp.ny;yi++){
		float y=(rp.recon_fov/rp.ny)*((float)yi-(rp.ny-1)/2.0f)+rp.y_origin;

		for (int k=0;k<n_half_turns;k++){
		    
		    // Geometric info about the ray (in-plane information)
		    float theta=tube_start+(2.0f*pi/cg.n_proj_turn)*i+k*pi;
		    float phat=x*sin(theta)-y*cos(theta);
		    float p_idx=phat/(cg.r_f*cg.fan_angle_increment/2.0f)+2.0f*cg.central_channel;

		    for (int zi=0;zi<mr->ri.n_slices_block;zi++){
			float z=zi*rp.coll_slicewidth+cg.z_rot/2.0f+cg.z_rot*tube_start/(2.0f*pi);

			// More geometric info about the ray (longitudinal information)
			float ray_pos=(cg.z_rot*(theta-asin(phat/cg.r_f))/(2.0f*pi));
			float lhat=sqrt(pow(cg.r_f,2.0f)-pow(phat,2.0f))-x*cos(theta)-y*sin(theta);
			float qhat=(z-ray_pos)/(lhat*tan(cg.theta_cone/2.0f));
			float q_idx=((qhat+1.0f)/2.0f)*(cg.n_rows-1.0f)+k*I*cg.n_rows;

			// Compute the voxel idx in output array. For output arrays, z is stored linearly in memory (stride 1),
			// then y (stride nz), then x (stride nz*ny)
			int vi=mr->ri.n_slices_block*mr->rp.ny*xi + mr->ri.n_slices_block*yi + zi;

			// Update voxel by interpolating from projection array and update weights
			s_t[vi]+=interp2(proj,dim,p_idx,q_idx)*W(qhat);
			h_t[vi]+=W(qhat);
		    }
		}
	    }
	}

	// Combine voxels with their weights and write back to final output array
	// In this array, x is stored linearly, y is stored with stride nx, z with stride nx*ny
	for (int xi=0;xi<mr->rp.nx;xi++){
	    for (int yi=0;yi<mr->rp.ny;yi++){
		for (int zi=0;zi<mr->ri.n_slices_block;zi++){
		    int out_idx=zi*mr->rp.nx*mr->rp.ny+yi*mr->rp.nx+xi;
		    int in_idx=mr->ri.n_slices_block*mr->rp.ny*xi + mr->ri.n_slices_block*yi + zi;

		    if (h_t[in_idx]!=0)
			h_output[out_idx]+=(1.0f/h_t[in_idx])*s_t[in_idx]*2*pi/cg.n_proj_turn;		
		}
	    }
	}

    }

    // Copy reconstructed sub-volume into the final full-size reconstructed volume
    long block_offset=(mr->ri.cb.block_idx-1)*mr->rp.nx*mr->rp.ny*mr->ri.n_slices_block;
    memcpy(&mr->ctd.image[block_offset],h_output,mr->rp.nx*mr->rp.ny*mr->ri.n_slices_block*sizeof(float));

    free(h_output);
    free(s_t);
    free(h_t);
    free(proj);
    
    return 0;
    
} // end backproject function
