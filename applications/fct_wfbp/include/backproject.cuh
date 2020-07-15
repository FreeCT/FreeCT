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

#define pi 3.1415368979f
#define K 1
#define I 16

texture<float,cudaTextureType2D,cudaReadModeElementType> tex_a;
texture<float,cudaTextureType2D,cudaReadModeElementType> tex_b;

__constant__ struct ct_geom d_cg;
__constant__ struct recon_params d_rp;

__device__ inline float W(float q){
    float out;
    float Q=0.6f;
    if (fabsf(q)<Q){
	out=1.0f;
    }
    else if ((fabsf(q)>=Q)&&(fabsf(q)<1.0f)){
	out=powf(cosf((pi/2.0f)*(fabsf(q)-Q)/(1.0f-Q)),2.0f);
    }
    else {
	out=0.0f;
    }
    return out;
}

__global__ void bp_a(float * output,int proj_idx,float tube_start,int n_half_turns){
    //Initialize sum voxels and weight voxels
    float s[K]={0};
    float s_t[K]={0};
    float h_t[K]={0};

    // Get xyz indices, convert to spatial coordinates
    int xi=threadIdx.x+blockIdx.x*blockDim.x;
    int yi=threadIdx.y+blockIdx.y*blockDim.y;
    int zi=K*(threadIdx.z+blockIdx.z*blockDim.z);

    float x;
    if (d_cg.table_direction==-1)
        x=(d_rp.recon_fov/d_rp.nx)*((float)xi-(d_rp.nx-1)/2.0f)+d_rp.x_origin;
    else
	x=(d_rp.recon_fov/d_rp.nx)*(-(float)xi+(d_rp.nx-1)/2.0f)+d_rp.x_origin;
    float y=(d_rp.recon_fov/d_rp.ny)*((float)yi-(d_rp.ny-1)/2.0f)+d_rp.y_origin;
    float z=zi*d_rp.coll_slicewidth+d_cg.z_rot/2.0f+d_cg.z_rot*tube_start/(2.0f*pi);
    
    for (int i=0;i<I;i++){	

	for (int ii=0;ii<K;ii++){ // zero out our holder arrays
	    s_t[ii]=0.0f;
	    h_t[ii]=0.0f;
	}
	
	for (int k=0;k<n_half_turns;k++){
	    float theta=tube_start+(2.0f*pi/d_cg.n_proj_turn)*(proj_idx+i)+k*pi;
	    float phat=x*sinf(theta)-y*cosf(theta);
	    float p_idx=phat/(d_cg.r_f*d_cg.fan_angle_increment/2.0f)+2.0f*d_cg.central_channel;
                
	    for (int j=0;j<K;j++){
		z+=j*d_rp.coll_slicewidth;
	    
		float ray_pos=(d_cg.z_rot*(theta-asinf(phat/d_cg.r_f))/(2.0f*pi));
		float lhat=sqrtf(d_cg.r_f*d_cg.r_f-phat*phat)-x*cosf(theta)-y*sinf(theta);
		float qhat=-d_cg.table_direction*(z-ray_pos)/(lhat*tanf(d_cg.theta_cone/2.0f));
		float q_idx=((qhat+1.0f)/2.0f)*(d_cg.n_rows-1.0f)+d_cg.n_rows*i+k*I*d_cg.n_rows;

		float interpolated_value=tex2D(tex_a,p_idx+0.5,q_idx+0.5)*W(qhat);

		if (!isnan(interpolated_value)){
		    s_t[j]+=interpolated_value;
		    h_t[j]+=W(qhat);
		}
		else{
		    h_t[j]=0.0f;
		}
		//s_t[j]+=tex2D(tex_a,p_idx+0.5,q_idx+0.5)*W(qhat);
		//h_t[j]+=W(qhat);

		//__syncthreads();
	    
	    }	    
	}
	
	for (int kk=0;kk<K;kk++){
	    if (h_t[kk]!=0){
		s[kk]+=(1.0f/h_t[kk])*s_t[kk];
	    }
	}
    }
    
    //__syncthreads();
    
    for (int k=0;k<K;k++){
	if (sqrt(x*x+y*y)<(d_cg.acq_fov/2.0f)){
          output[d_rp.nx*d_rp.ny*(zi+k)+d_rp.nx*xi+yi]+=s[k]*2*pi/d_cg.n_proj_turn;//(1.0f/h[k])*s[k];        
	}
	else{
	    output[d_rp.nx*d_rp.ny*(zi+k)+d_rp.nx*xi+yi]=0.0f;
	}
    }
    
}

__global__ void bp_b(float * output,int proj_idx,float tube_start,int n_half_turns){
    //Initialize sum voxels and weight voxels
    float s[K]={0};
    float s_t[K]={0};
    float h_t[K]={0};
    
    // Get xyz indices, convert to spatial coordinates
    int xi=threadIdx.x+blockIdx.x*blockDim.x;
    int yi=threadIdx.y+blockIdx.y*blockDim.y;
    int zi=K*(threadIdx.z+blockIdx.z*blockDim.z);

    float x;
    if (d_cg.table_direction==-1)
	x=(d_rp.recon_fov/d_rp.nx)*((float)xi-(d_rp.nx-1)/2.0f)+d_rp.x_origin;
    else
	x=(d_rp.recon_fov/d_rp.nx)*(-(float)xi+(d_rp.nx-1)/2.0f)+d_rp.x_origin;
    float y=(d_rp.recon_fov/d_rp.ny)*((float)yi-(d_rp.ny-1)/2.0f)+d_rp.y_origin;
    float z=zi*d_rp.coll_slicewidth+d_cg.z_rot/2.0f+d_cg.z_rot*tube_start/(2.0f*pi);
    
    for (int i=0;i<I;i++){
	
	for (int ii=0;ii<K;ii++){ // zero out our holder arrays
	    s_t[ii]=0.0f;
	    h_t[ii]=0.0f;
	}
	
	for (int k=0;k<n_half_turns;k++){
	    float theta=tube_start+(2.0f*pi/d_cg.n_proj_turn)*(proj_idx+i)+k*pi;
	    float phat=x*sinf(theta)-y*cosf(theta);
	    float p_idx=phat/(d_cg.r_f*d_cg.fan_angle_increment/2.0f)+2.0f*d_cg.central_channel;
	    
	    for (int j=0;j<K;j++){
		z+=j*d_rp.coll_slicewidth;
	    
		float ray_pos=(d_cg.z_rot*(theta-asin(phat/d_cg.r_f))/(2.0f*pi));
		float lhat=sqrt(d_cg.r_f*d_cg.r_f-phat*phat)-x*cosf(theta)-y*sinf(theta);
		float qhat=-d_cg.table_direction*(z-ray_pos)/(lhat*tanf(d_cg.theta_cone/2.0f));
		float q_idx=((qhat+1.0f)/2.0f)*(d_cg.n_rows-1.0f)+d_cg.n_rows*i+k*I*d_cg.n_rows;

		float interpolated_value=tex2D(tex_b,p_idx+0.5,q_idx+0.5)*W(qhat);

		if (!isnan(interpolated_value)){
		    s_t[j]+=interpolated_value;
		    h_t[j]+=W(qhat);
		}
		else{
		    h_t[j]=0.0f;
		}
		
		//s_t[j]+=tex2D(tex_b,p_idx+0.5,q_idx+0.5)*W(qhat);
		//h_t[j]+=W(qhat);
		    
		//__syncthreads();

	    }
	}
	
	for (int kk=0;kk<K;kk++){
	    if (h_t[kk]!=0){
		//if ((h_t[kk]!=0)&&(!isnan(s[kk]))){		
		s[kk]+=(1.0f/h_t[kk])*s_t[kk];
	    }
	}
    }

    //__syncthreads();
    
    for (int k=0;k<K;k++){
	if (sqrt(x*x+y*y)<(d_cg.acq_fov/2.0f)){
          output[d_rp.nx*d_rp.ny*(zi+k)+d_rp.nx*xi+yi]+=s[k]*2*pi/d_cg.n_proj_turn;//(1.0f/h[k])*s[k];
	}
	else{
	    output[d_rp.nx*d_rp.ny*(zi+k)+d_rp.nx*xi+yi]=0.0f;//(1.0f/h[k])*s[k];
	}
    }    
}


