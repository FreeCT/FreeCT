/* FreeCT_ICD is MBIR CT reconstruction Software */
/* Copyright (C) 2018  John Hoffman, Stefano Young, Frederic Noo */

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
/* jmhoffman@mednet.ucla.edu with "FreeCT_ICD" in the subject line*/

#include <iostream>
#include <cstring>
#include <cmath>

#include "spinner.h"
#include "recon_structs.h"
#include "interp.h"

#define PI 3.141592653589793238

// From Xu, Tsui 2012:
//   (*) the slice thickness dz = l*h/N_{2pi}, l > 0 is an integer; 
//   (*) the rotation angle between adjacent slices is l*2*pi/N_{2pi},
//
//    h is table feed per rotation, N_{2pi} is number of views per rotation

void rotate_slice(float * in_slice, float * out_slice, double angle,const struct recon_params * rp);

inline bool in(double x, double min, double max){
    bool tf=false;
    
    if (x>=min && x<=max)
        tf=true;

    return tf;
}

void rotate_slices_fixed2rotating(const struct recon_params * rp, struct ct_data * data){    
    // Rotate slices from fixed cartesian grid into rotating coordinate frame
    // e.g. Initial WFBP slices into ICD rotating frame

    std::cout << "Rotating WFBP slices to match ICD w/ rotating grid" << std::endl;

    // Allocate a temporary buffer to hold slice during rotation
    size_t n_pixels=rp->num_voxels_x*rp->num_voxels_y;
    float * tmp_slice=new float[n_pixels];

    init_spinner();
    
    for (int i=0; i<rp->num_voxels_z; i++){

        update_spinner((size_t)i,rp->num_voxels_z);
        
        // Determine the amount we need to rotate (based on slice index)
        double rotation_angle=-(data->tube_angles[data->slice_indices[i]])*2.0*PI/360.0;
        rotation_angle=-rotation_angle;
        
        // Copy current slice to temporary buffer
        std::memcpy((char*)tmp_slice,(char*)&data->recon_volume[i*n_pixels],n_pixels*sizeof(float));
        
        // Rotate slice back into our array
        rotate_slice(tmp_slice,&data->recon_volume[i*n_pixels],rotation_angle,rp);
    }

    destroy_spinner();
    
    // Free allocated memory
    delete[] tmp_slice;
}

void rotate_slices_rotating2fixed(const struct recon_params * rp, struct ct_data * data){
    // Return rotating slices to fixed cartesian grid
    // e.g. Reconstructed ICD slices back for final save to disk

    std::cout << "Rotating ICD slices into fixed cartesian grid" << std::endl;

    // Allocate a temporary buffer to hold slice during rotation
    size_t n_pixels=rp->num_voxels_x*rp->num_voxels_y;
    float * tmp_slice=new float[n_pixels];

    init_spinner();

    for (int i=0; i<rp->num_voxels_z; i++){

        update_spinner((size_t)i,rp->num_voxels_z);

        //std::cout << "Slice Index: " << i << ": "<<   data->slice_indices[i] << std::endl;

        // Determine the amount we need to rotate (based on slice index)
        double rotation_angle=-(data->tube_angles[data->slice_indices[i]])*2.0*PI/360.0;
        
        // Copy current slice to temporary buffer
        std::memcpy((char*)tmp_slice,(char*)&data->recon_volume[i*n_pixels],n_pixels*sizeof(float));
        
        // Rotate slice back into our array
        rotate_slice(tmp_slice,&data->recon_volume[i*n_pixels],rotation_angle,rp);
    }

    destroy_spinner();
    
    // Free allocated memory
    delete[] tmp_slice;
    
}

void rotate_slice(float * slice, float * out_slice,double angle,const struct recon_params * rp){

    // Rotates based on COORDINATES ONLY (i.e. does not transform to spatial coordinates)
    // Assumes that the center for reconstructed volume is the center of gantry

    // "angle" has units of radians
    
    // Step through voxels in the output matrix
    for (int x_idx=0; x_idx<rp->num_voxels_x; x_idx++){
        for (int y_idx=0; y_idx<rp->num_voxels_y; y_idx++){
            
            // Determine the output index
            size_t out_idx=x_idx+y_idx*rp->num_voxels_y;
            
            // Determine coordinates of current_voxel in output matrix rel. to center
            double x = x_idx - rp->center_voxel_x;
            double y = y_idx - rp->center_voxel_y;

            // Apply rotation matrix
            // [[x']  =  [[ cos(th) , -sin(th)]  * [[x]
            //  [y']] =   [ sin(th) ,  cos(th)]] *  [y]]
            //
            // =>
            //   x'=cos(th)*x-sin(th)*y
            //   y'=sin(th)*x+cos(th)*y

            float x_prime=cos(angle)*x-sin(angle)*y;
            float y_prime=sin(angle)*x+cos(angle)*y;


            // Interpolate from original image and store in the output image
            struct array_dims dims;
            dims.idx1=rp->nx;
            dims.idx2=rp->ny;
            
            out_slice[out_idx]=interp2(slice,dims,x_prime+rp->center_voxel_x,y_prime+rp->center_voxel_y);
        }        
    }    
}
