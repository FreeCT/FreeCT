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

/* This file was automatically generated using a code-generation
   script.  Changes made directly to the .c file will likely be
   overwritten.  If you have suggested changes to the script we may
   implement them in the templates.  Please contact
   freect.project@gmail.com.*/

#pragma once
#ifndef parse_config_h
#define parse_config_h

#include <stdlib.h>
#include <stdio.h>
#include <string>

#define COMMENT_DELIM '%'

struct recon_params {
        char raw_data_dir[4096];
    char raw_data_file[255];
    char output_dir[4096];
    char output_file[255];
    int n_rows;
    float coll_slicewidth;
    float start_pos;
    float end_pos;
    float pitch_value;
    float slice_thickness;
    float acq_fov;
    float recon_fov;
    int recon_kernel;
    int n_readings;
    float x_origin;
    float y_origin;
    int z_ffs;
    int phi_ffs;
    char scanner[4351];
    int file_type;
    int file_subtype;
    int raw_data_offset;
    unsigned int nx;
    unsigned int ny;
    float tube_start_angle;
    float adaptive_filtration_s;
    int n_slices;
    char table_dir_str[1024];
    int table_dir;

};

static inline void empty_config(const char * filepath){
    char fullpath[4096+255]={0};
    strcpy(fullpath,filepath);

    FILE * fid = fopen(fullpath,"w");

    fprintf(fid,

    "RawDataDir:\n"
"RawDataFile:\n"
"OutputDir:\n"
"OutputFile:\n"
"Nrows:\n"
"CollSlicewidth:\n"
"StartPos:\n"
"EndPos:\n"
"TableFeed:\n"
"PitchValue:\n"
"SliceThickness:\n"
"AcqFOV:\n"
"ReconFOV:\n"
"ReconKernel:\n"
"Readings:\n"
"Xorigin:\n"
"Yorigin:\n"
"Zffs:\n"
"Phiffs:\n"
"Scanner:\n"
"FileType:\n"
"FileSubType:\n"
"RawOffset:\n"
"Nx:\n"
"Ny:\n"
"TubeStartAngle:\n"
"AdaptiveFiltration:\n"
"NSlices:\n"
"TableDir:\n"
"TableDirInt:\n"


    );

    fclose(fid);
    exit(0);
};

//void parse_config(char * config_file, struct recon_params * structure);
void parse_config(std::string config_file, struct recon_params * structure);


#endif