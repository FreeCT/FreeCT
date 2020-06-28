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

#ifndef interp_h
#define interp_h

#include <math.h>

struct array_dims{
    int idx1;
    int idx2;
    int idx3;
};

float interp1(float * array, float idx1);
float interp2(float * array, struct array_dims dim, float idx1,float idx2);
float interp3(float * array, struct array_dims dim, float idx1,float idx2,float idx3);

float interp1(float * array, float idx1){
    float w=idx1-floor(idx1);
    return array[(int)idx1]*(1.0f-w)+array[(int)idx1+1]*w;
}

float interp2(float * array, struct array_dims dim, float idx1,float idx2){
    //Assumes idx1 is stored linearly in memory with stride 1
    //idx2 is stored with stride dim_idx 1
    //
    //     dim1------->
    //dim2  0   0   0   0   0   0   0   0
    //|     1   1   1   1   1   1   1   1
    //|     2   2   2   2   2   2   2   2 
    //|     .           .
    //|     .               .
    //|     .                   .

    float val;

    // Edge mode "Border"
    // i.e. if user requests data outside of array, return 0; 
    if ( idx1>(dim.idx1-2) || (idx1<0) || idx2>(dim.idx2-2) || (idx2<0) )
	val=0.0f;
    else {
    float v=idx1-floor(idx1);
    float w=idx2-floor(idx2);

    int a_idx1=(int)floor(idx2)*dim.idx1+(int)floor(idx1);
    int a_idx2=(int)floor(idx2)*dim.idx1+(int)floor(idx1)+1;
    int a_idx3=((int)floor(idx2)+1)*dim.idx1+(int)floor(idx1);
    int a_idx4=((int)floor(idx2)+1)*dim.idx1+(int)floor(idx1)+1;

    val=   array[a_idx1]  *  (1.0f-v)  *(1.0f-w)
	+  array[a_idx2]  *      v     *(1.0f-w)
	+  array[a_idx3]  *  (1.0f-v)  *    w
	+  array[a_idx4]  *      v     *    w;

    }
    
    return val;
}

float interp3(float * array, struct array_dims dim, float idx1,float idx2,float idx3){
    // dim.idx1 and idx1 are stored linearly in memory
    // dim.idx2 and idx2 are stored in memory with stride dim.idx1
    // dim.idx3 and idx3 are stored in memory with stride dim.idx1*dim.idx2
    
    // Clamping
    if (idx1>(dim.idx1-1))
	idx1=dim.idx1-1.0f;
    if (idx1<0)
	idx1=0.0f;
    if (idx2>(dim.idx2-1))
	idx2=dim.idx2-1.0f;
    if (idx2<0)
	idx2=0.0f;
    if (idx3>(dim.idx3-1))
	idx3=dim.idx3-1.0f;
    if (idx3<0)
	idx3=0.0f;

    // Find weights
    float u=idx1-floor(idx1);
    float v=idx2-floor(idx2);
    float w=idx3-floor(idx3);

    // Find linear indices for interpolation points
    int a_idx1=(int)floor(idx3)*dim.idx2*dim.idx1     + (int)floor(idx2)*dim.idx1     + (int)floor(idx1);
    int a_idx2=(int)floor(idx3)*dim.idx2*dim.idx1     + (int)floor(idx2)*dim.idx1     + (int)floor(idx1) + 1;
    int a_idx3=(int)floor(idx3)*dim.idx2*dim.idx1     + ((int)floor(idx2)+1)*dim.idx1 + (int)floor(idx1);
    int a_idx4=(int)floor(idx3)*dim.idx2*dim.idx1     + ((int)floor(idx2)+1)*dim.idx1 + (int)floor(idx1) + 1;    
    int a_idx5=((int)floor(idx3)+1)*dim.idx2*dim.idx1 + (int)floor(idx2)*dim.idx1     + (int)floor(idx1);
    int a_idx6=((int)floor(idx3)+1)*dim.idx2*dim.idx1 + (int)floor(idx2)*dim.idx1     + (int)floor(idx1) + 1;
    int a_idx7=((int)floor(idx3)+1)*dim.idx2*dim.idx1 + ((int)floor(idx2)+1)*dim.idx1 + (int)floor(idx1);
    int a_idx8=((int)floor(idx3)+1)*dim.idx2*dim.idx1 + ((int)floor(idx2)+1)*dim.idx1 + (int)floor(idx1) + 1;    

    // Clamp any indices that would be over the edge of our array
    if (idx1==dim.idx1-1.0f){
	a_idx2=a_idx1;
	a_idx4=a_idx3;
	a_idx6=a_idx5;
	a_idx8=a_idx7;
    }
    if (idx2==dim.idx2-1.0f){
	a_idx3=a_idx1;
	a_idx4=a_idx2;
	a_idx7=a_idx5;
	a_idx8=a_idx6;

    }
    if (idx3==dim.idx3-1.0f){
	a_idx5=a_idx1;
	a_idx6=a_idx2;
	a_idx7=a_idx3;
	a_idx8=a_idx4;
    }
    
    //Return the interpolation
    return array[a_idx1]  *  (1-u) * (1-v) * (1-w) + 
	   array[a_idx2]  *    u   * (1-v) * (1-w) +
	   array[a_idx3]  *  (1-u) *   v   * (1-w) +
	   array[a_idx4]  *    u   *   v   * (1-w) +
	   array[a_idx5]  *  (1-u) * (1-v) *   w   +
	   array[a_idx6]  *    u   * (1-v) *   w   +
	   array[a_idx7]  *  (1-u) *   v   *   w   +
           array[a_idx8]  *    u   *   v   *   w   ;
}


#endif
