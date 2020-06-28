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

#ifndef recon_structs_h
#define recon_structs_h

#include <parse_config.h>

struct block_info{
    int block_idx;
    float block_slice_start;
    float block_slice_end;
    int idx_block_slice_start;
    int idx_block_slice_end;
};

struct recon_info{
    size_t n_ffs;
    float data_begin_pos;
    float data_end_pos;
    float allowed_begin;
    float allowed_end;
    size_t n_slices_requested;
    size_t n_slices_recon;
    size_t n_slices_block;
    size_t n_blocks;
    size_t idx_slice_start;
    size_t idx_slice_end; 
    float recon_start_pos;
    float recon_end_pos;
    size_t idx_pull_start;
    size_t idx_pull_end;
    size_t n_proj_pull;
    struct block_info cb;
};
    
struct ct_data{
    float * raw;
    float * rebin;
    float * image;
    float * d_raw;
    float * d_rebin;
    float * d_image;
    float * d_final_image_stack;
    float * final_image_stack;
};
    
struct ct_geom{
    size_t n_proj_turn;
    size_t n_proj_ffs;
    size_t n_channels;
    size_t n_channels_oversampled;
    size_t n_rows;
    size_t n_rows_raw;
    float r_f;
    float z_rot;
    float theta_cone;
    float fan_angle_increment;
    float src_to_det;
    float anode_angle;
    float central_channel;
    float acq_fov;
    size_t projection_offset;
    size_t add_projections;
    size_t add_projections_ffs;
    int reverse_row_interleave;
    int reverse_channel_interleave;

    // -1 table positions decreasing (SciDirTableIn);
    //  1 table positions increasing (SciDirTableOut);
    int table_direction;
};

struct flags{
    int testing;
    int verbose;
    int no_gpu;
    int set_device;
    int device_number;
    int timing;
    int benchmark;
};
    
struct recon_metadata {
    char home_dir[4096];
    char cwd[4096];
    char install_dir[4096];
    struct flags flags;
    struct recon_params rp;
    struct recon_info ri;
    struct ct_geom cg;
    struct ct_data ctd;
    
    float * tube_angles;
    double * table_positions;
};

#endif
