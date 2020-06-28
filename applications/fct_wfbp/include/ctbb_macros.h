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

#ifndef ctbb_macros_h
#define ctbb_macros_h

/* GPU Error Handling */
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
    if (code != cudaSuccess){
	fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
	if (abort) exit(code);
    }
}

/* Debugging */
#define FLOAT_DEBUG(array_ptr,numel,file) float_debug(array_ptr,numel,file);
inline void float_debug(float * array, size_t numel, const char * filename){
    FILE * fid=fopen(filename,"w");
    fwrite(array,sizeof(float),numel,fid);
    fclose(fid);
}

/* Timing */
// Mini macros (you probably don't need to use these directly
#define __TIMING_VARIABLE__ milli
#define TIMER_INIT() cudaEvent_t start,stop; float __TIMING_VARIABLE__;
#define TIMER_START() {cudaEventCreate(&start);cudaEventCreate(&stop);cudaEventRecord(start);}
#define TIMER_END() {cudaEventRecord(stop);			\
	cudaEventSynchronize(stop);				\
	cudaEventElapsedTime(&__TIMING_VARIABLE__,start,stop);	\
	cudaEventDestroy(start);				\
	cudaEventDestroy(stop);}
#define TIMER_PRINT(flag,message) print_time(flag,message,__TIMING_VARIABLE__);
inline void print_time(int flag, const char * message,float time){
    if (flag)
	printf("%0.2f ms for %s\n",time,message);
}

#define __MASTER_TIMING_VARIABLE__ m_milli
#define TIMER_MASTER_INIT() cudaEvent_t m_start,m_stop; float __MASTER_TIMING_VARIABLE__ = 0.0f;
#define TIMER_MASTER_START() cudaEventCreate(&m_start);cudaEventCreate(&m_stop);cudaEventRecord(m_start);
#define TIMER_MASTER_END() cudaEventRecord(m_stop);			\
    cudaEventSynchronize(m_stop);					\
    cudaEventElapsedTime(&__MASTER_TIMING_VARIABLE__,m_start,m_stop);	\
    cudaEventDestroy(m_start);						\
    cudaEventDestroy(m_stop);
#define TIMER_MASTER_PRINT(flag,message) print_time(flag,message,__MASTER_TIMING_VARIABLE__);


// This is what we use to time a function call
#define TIME_EXEC(function_call,flag,message) TIMER_START(); function_call; TIMER_END(); TIMER_PRINT(flag,message);

#endif
