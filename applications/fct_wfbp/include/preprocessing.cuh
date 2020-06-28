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

#define F_MAX 0.03f
#define S 1.0f
#define pi 3.1415368979f

#define __FILTER_3D__ //__FILTER_3D__

__constant__ struct ct_geom d_cg;
__constant__ struct recon_info d_ri;

__device__ float rhs(float * array,int start_idx, int max_idx,int numel_range,float T){


    /* NOTE:!!!!!!!! When attentuations are high, this truncating
       seems to have the effect of reducing overall thresholds
       drastically, which will likely cause oversmoothing.  To combat
       this, we may want to clamp any value that can't get a good,
       full range to some value.  Currently however, we've left it
       as is.  There could also be an error.*/

    int sum=0;

    // Handle edge cases (i.e.)
    int diff=start_idx+(numel_range-1);
    if (diff<(numel_range-1)){
	numel_range+=start_idx;
	start_idx=0;
    } // spilling over lhs of array (start idx<0)

    if (diff>max_idx){ 
	numel_range-=(diff-max_idx);
    } // spilling over rhs side of array

    
    for (int i=0; i<numel_range; i++){
	if (array[start_idx+i]>T)
	    sum+=1;
    }

    return (1.0f/pi)*(float)sum*(1.0f/(float)d_cg.n_proj_ffs);
}

__device__ float triangle_3d(int N, int i, int j, int k){
    int cN2=(int) ceilf((float)N/2.0f);
    return (float)((cN2-abs(i-cN2))+(cN2-abs(j-cN2))+(cN2-abs(k-cN2))+1); 
}

__global__ void extract_sup(float * raw,float * sup){

    size_t proj_idx=threadIdx.x+blockIdx.x*blockDim.x;
    size_t proj_size=d_cg.n_channels*d_cg.n_rows_raw;
    size_t row_size=d_cg.n_channels;

    float max_val=0.0f;
    
    for (int i=0; i<d_cg.n_channels; i++){
	for (int j=0; j<d_cg.n_rows_raw; j++){
	    max_val=max(max_val,raw[proj_size*proj_idx+j*row_size+i]);
	}
    }

    sup[proj_idx]=max_val;
    
}

__global__ void smooth_sup(float * sup_raw,float * sup_smooth){

    int proj_idx=threadIdx.x+blockIdx.x*blockDim.x;
    
    // 18 degree value pulled from Kachelreiss and Kalendar 2001
    int n_avg=ceilf((18.0f/360.0f)*(float)d_cg.n_proj_ffs);
    n_avg=n_avg+(1-n_avg%2);
    
    float running_avg=0.0;
    int offset=floorf(n_avg/2);
    int max_idx=d_ri.n_proj_pull-1;
    for (int i=0; i<n_avg;i++){

	int raw_idx=proj_idx-offset+i;

	//Clamp
	if (raw_idx<0)
	    raw_idx=0;
	else if (raw_idx>max_idx)
	    raw_idx=max_idx;

	running_avg+=(1.0f/(float)n_avg)*sup_raw[raw_idx];
    }
    
    sup_smooth[proj_idx]=running_avg;

}

__global__ void eccentricity(float * sup_smooth, float * ecc,float * max_array,float * min_array){
    int proj_idx=threadIdx.x+blockIdx.x*blockDim.x;
    
    // 180 degree value pulled from Kachelreiss and Kalendar 2001
    int n_range=ceilf((180.0f/360.0f)*(float)d_cg.n_proj_ffs);
    n_range=n_range+(1-n_range%2);
    
    int offset=floorf(n_range/2);
    int max_idx=d_ri.n_proj_pull-1;

    float max_val=0.0f;
    float min_val=10000.0f;
    //Find max and min over 180 degrees
    for (int i=0; i<n_range;i++){

	int raw_idx=proj_idx-offset+i;

	//Clamp
	if (raw_idx<0)
	    raw_idx=0;
	else if (raw_idx>max_idx)
	    raw_idx=max_idx;

	float v=sup_smooth[raw_idx];
	
	if (v>max_val)
	    max_val=v;

	if (v<min_val)
	    min_val=v;
	
    }

    // Calculate raw eccentricity
    max_array[proj_idx]=max_val;
    min_array[proj_idx]=min_val;
    
    float e=1-min_val/max_val;

    //Rescale range [0.3 0.5] to [0 1] (values from Kachelreiss Kalendar 2001)
    e=(e-0.3f)/(0.5f-0.3f);

    if (e>1)
	ecc[proj_idx]=1.0f;
    else if (e<0)
	ecc[proj_idx]=0.0f;
    else
	ecc[proj_idx]=e;
}

__global__ void eccentricity_rescale(float * ecc, float * ecc_rescale){}

__global__ void find_thresholds(float strength,float * ecc_rescale,float * sup_smooth, float * max_array, float * min_array,float * thresholds){
    int proj_idx=threadIdx.x+blockIdx.x*blockDim.x;
    
    // 180 degree value pulled from Kachelreiss and Kalendar 2001
    int n_range=ceilf((180.0f/360.0f)*(float)d_cg.n_proj_ffs);
    n_range=n_range+(1-n_range%2);
    
    int offset=floorf(n_range/2);
    int max_idx=d_ri.n_proj_pull-1;
    int max_iter=20;

    // Start values
    float lower_bound=min_array[proj_idx];
    float upper_bound=max_array[proj_idx];
    float T=(max_array[proj_idx]+min_array[proj_idx])/2.0f;
    float lhs_alpha=F_MAX*strength*ecc_rescale[proj_idx];
    float rhs_alpha=rhs(sup_smooth,proj_idx-offset,max_idx,n_range,T);
    
    // Binary search over 
    int iter=0;
    while ((abs(lhs_alpha-rhs_alpha)/lhs_alpha)>0.01){

	if (iter>=max_iter)
	    break;
	else
	    iter+=1;

	if (rhs_alpha>lhs_alpha)
	    lower_bound=T;
	else if (rhs_alpha<lhs_alpha)
	    upper_bound=T;
	else
	    break;

	T=(lower_bound+upper_bound)/2.0f;
	rhs_alpha=rhs(sup_smooth,proj_idx-offset,max_idx,n_range,T);
    }

    thresholds[proj_idx]=T;
}

__global__ void filter_projections(float * raw,float * thresholds, float * filtered_raw){

    size_t proj_idx=threadIdx.x+blockIdx.x*blockDim.x;
    int n_proj=d_ri.n_proj_pull;
    int n_channels=d_cg.n_channels;
    int n_rows_raw=d_cg.n_rows_raw;
    size_t proj_size=n_channels*n_rows_raw;
    size_t row_size=n_channels;

    float T=thresholds[proj_idx];

#ifndef __FILTER_3D__
    float triangle_2d[25]={1,2,3,2,1, 2,3,4,3,2, 3,4,5,4,3, 2,3,4,3,2, 1,2,3,2,1};
#endif

    
    for (int i=0; i<n_channels; i++){
	for (int j=0; j<n_rows_raw; j++){

	    size_t idx=proj_size*proj_idx+j*row_size+i;

	    if (raw[idx]>T){
		
#ifdef __FILTER_3D__
		// 3D Filter
		float norm=0.0f;
		float filtered_val=0.0f;
		for (int ii=0;ii<5;ii++){//loop over filter channels
		    for (int jj=0;jj<5;jj++){//loop over filter rows
			for (int kk=0;kk<5;kk++){//loop over filter projections

			    int f_row_idx=j+jj-2; // proj_row_idx+filter_row_idx - offset
			    int f_channel_idx=i+ii-2;
			    int f_proj_idx=proj_idx+kk-2;
			    
			    if ((f_row_idx<0||f_channel_idx<0||f_proj_idx<0)||
				(f_row_idx>=n_rows_raw||f_channel_idx>=n_channels||f_proj_idx>=n_proj)){
				continue;
			    }
			    else{
				size_t idx_tmp=proj_size*f_proj_idx+row_size*f_row_idx+f_channel_idx;
				float triangle_val=triangle_3d(5,ii,jj,kk);
				norm+=triangle_val;
				filtered_val+=raw[idx_tmp]*triangle_val;
			    }			    
			}
		    }
		}
		
#else
		// 2D filter
		float norm=0.0f;
		float filtered_val=0.0f;
		for (int ii=0;ii<5;ii++){//loop over filter channels
		    for (int jj=0;jj<5;jj++){//loop over filter rows
			int f_row_idx=j+jj-2; // proj_row_idx+filter_row_idx - offset
			int f_channel_idx=i+ii-2;
			
			if ((f_row_idx<0||f_channel_idx<0)||(f_row_idx>=n_rows_raw||f_channel_idx>=n_channels)){
			    continue;
			}
			else{
			    size_t idx_tmp=proj_size*proj_idx+f_row_idx*row_size+f_channel_idx;
			    norm+=triangle_2d[jj*5+ii];
			    filtered_val+=raw[idx_tmp]*triangle_2d[jj*5+ii];
			}
		    }
		}
#endif
		filtered_raw[idx]=(1.0f/norm)*filtered_val;
	    }
	    else{
		filtered_raw[idx]=raw[idx];
	    }
	}
    }

} 
