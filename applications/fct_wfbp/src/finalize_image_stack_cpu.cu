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
#include <finalize_image_stack_cpu.h>

int finalize_image_stack_cpu(struct recon_metadata * mr){

    struct recon_params rp=mr->rp;
    struct recon_info ri=mr->ri;

    float * temp_out=(float*)calloc(rp.nx*rp.ny*ri.n_slices_recon,sizeof(float));

    int recon_direction=fabs(rp.end_pos-rp.start_pos)/(rp.end_pos-rp.start_pos);
    if (recon_direction!=1&&recon_direction!=-1) // user request one slice (end_pos==start_pos)
	recon_direction=1;
    
    int table_direction=(mr->table_positions[1000]-mr->table_positions[0])/abs(mr->table_positions[1000]-mr->table_positions[0]);
    if (table_direction!=1&&table_direction!=-1)
	printf("Axial scans are currently unsupported, or a different error has occurred\n");
    
    // Check for a reversed stack of images and flip, otherwise just copy
    if (recon_direction!=table_direction){
	for (int b=0;b<ri.n_blocks;b++){
	    for (int z=0;z<ri.n_slices_block;z++){
		for (int x=0;x<rp.nx;x++){
		    for (int y=0;y<rp.ny;y++){
			long block_offset=b*rp.nx*rp.ny*ri.n_slices_block;
			temp_out[z*rp.nx*rp.ny+y*rp.nx+x+block_offset]=mr->ctd.image[((ri.n_slices_block-1)-z)*rp.nx*rp.ny+y*rp.nx+x+block_offset];
		    }
		}
	    }
	}
    }
    else{
	for (int z=0;z<ri.n_slices_recon;z++){
	    for (int x=0;x<rp.nx;x++){
		for (int y=0;y<rp.ny;y++){
		    temp_out[z*rp.nx*rp.ny+y*rp.nx+x]=mr->ctd.image[z*rp.nx*rp.ny+y*rp.nx+x];
		}
	    }
	}	
    }

    // Once we have straightened our image stack out, we need to adjust slice thickness
    // to match what the user requested.
    // We use a triangle average with the FWHM equal to the requested slice thickness

    int n_raw_images=ri.n_slices_block*ri.n_blocks;
    int n_slices_final=floor(fabs(rp.end_pos-rp.start_pos)/rp.slice_thickness)+1;
    mr->ctd.final_image_stack=(float*)calloc(rp.nx*rp.ny*n_slices_final,sizeof(float));

    float * recon_locations;
    recon_locations=(float*)calloc(n_slices_final,sizeof(float));
    for (int i=0;i<n_slices_final;i++){
	recon_locations[i]=rp.start_pos+recon_direction*i*rp.slice_thickness;
    }

    float * raw_recon_locations;
    raw_recon_locations=(float*)calloc(n_raw_images,sizeof(float));
    for (int i=0;i<n_raw_images;i++){
	raw_recon_locations[i]=ri.recon_start_pos+recon_direction*i*rp.coll_slicewidth;//(rp.start_pos-recon_direction*rp.slice_thickness)+recon_direction*i*rp.coll_slicewidth;
    }

    float * weights;
    weights=(float*)calloc(n_raw_images,sizeof(float));

    // Loop over slices
    for (int k=0;k<n_slices_final;k++){
	float slice_location=recon_locations[k];
	// Calculate all of the weights for the unaveraged slices
	float sum_weights=0;
	for (int step=0;step<ri.n_slices_block*ri.n_blocks;step++){
	    weights[step]=fmax(0.0f,1.0f-fabs(raw_recon_locations[step]-slice_location)/rp.slice_thickness);
	    sum_weights+=weights[step];
	}
	
	// Loop over pixels in slice k
	for (int i=0;i<rp.nx;i++){
	    for (int j=0;j<rp.ny;j++){
		// Carry out the averaging
		int out_idx=k*rp.nx*rp.ny+j*rp.nx+i;
		for (int raw_slice=0;raw_slice<n_raw_images;raw_slice++){		    
		    if (weights[raw_slice]!=0){
			int raw_idx=raw_slice*rp.nx*rp.ny+j*rp.nx+i;
			mr->ctd.final_image_stack[out_idx]+=(weights[raw_slice]/sum_weights)*temp_out[raw_idx];
		    }
		}
	    }
	}
    }

    // Free stuff allocated inside cleanup
    free(temp_out);
    free(recon_locations);
    free(raw_recon_locations);
    free(weights);

    return 0;
}
