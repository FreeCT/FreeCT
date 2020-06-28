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
#include <fstream>
#include <cstring>
#include <string>
#include <cstdlib>

#include "recon_structs.h"
#include "initialize_recon_volume.h"
#include "rotate_slices.h"

//#define __WFBP_RECON__ "fct_wfbp" //Must be either "ctbb_recon" or "fct_wfbp"
#define __WFBP_RECON__ "ctbb_recon" //Must be either "ctbb_recon" or "fct_wfbp"

void build_wfbp_parameter_file(const struct recon_params * rp);

inline bool exists(const std::string& name){
    std::ifstream f(name.c_str());
    return f.good();
}

inline std::string split_path(std::string filepath){
    size_t idx=filepath.find_last_of('/');
    return filepath.substr(0,idx);    
}

inline std::string split_filename(std::string filepath){
    size_t idx=filepath.find_last_of('/');
    return filepath.substr(idx+1);    
}

void initialize_recon_volume(const struct recon_params * rp, struct ct_data * data){

    /* Load an existing reconstruction from a file */
    // No validation is performed on the file to decided if it's
    // a reconstruction, or from the same dataset, or rotated properly.    
    // User is responsible for guaranteeing that all of these things are true
    if (rp->initial_recon_path.compare("")!=0){        
        std::cout << "Loading initial reconstruction from file" << std::endl;

        if (!exists(rp->initial_recon_path)){
            std::cout << "Specified initial reconstruction file not found. Exiting" << std::endl;
            exit(1);
        }

        std::ifstream infile(rp->initial_recon_path,std::ios::binary);

        infile.read((char*)data->recon_volume,rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z*sizeof(float));

        if (!infile){
            std::cout << "File was not read properly." << std::endl;
            std::cout << "Could not initialize reconstruction volume from given file. Exiting" << std::endl;
            exit(1);
        }

        infile.close();
    }

    /* Initialize the recontruction volume with a WFBP reconstruction */
    // This will handle the initial reconstruction across the requested range
    // as well as ensure that the slices are rotated properly to maximize the
    // utility of the WFBP recon to accelerate convergence.
    else if (rp->wfbp_initialize!=0){

        // Check that FreeCT/CTBangBang is installed on the machine
        // (To do this, we make a system call and check exit
        // status. FreeCT/CTBangBang will return a usage statement and
        // exit code 0 if run without any arguments. This is not
        // especially robust and likely has edge cases where it may
        // fail, however if all of the software is installed using
        // "defaults" this is a fine approach.)

        int status=system(__WFBP_RECON__);

        if (status){
            std::cout << "WFBP reconstruction software, " <<  __WFBP_RECON__ << ", not found." << std::endl;
            std::cout << "WFBP recon software is available at: http://cvib.ucla.edu/freect" << std::endl;
            std::cout << "To start from blank volume, comment out line \"wfbp_initialize: ...\" in parameter file. (Much longer execution time)." << std::endl;
            std::cout << "\n" << std::endl;
            std::cout << "Unable to initialize from WFBP recon." << std::endl;
            std::cout << "Exiting." << std::endl;
            exit(1);
        }

        // Build the input parameter file from existing ICD parameters
        build_wfbp_parameter_file(rp);
        
        // Call the FreeCT program on the built parameter file
        std::string program=__WFBP_RECON__;
        std::string options=" -v --timing ";
        std::string parameter_file=rp->output_dir + "/init_wfbp.prm";

        std::string recon_command = program + options + parameter_file;
        
        status=system(recon_command.c_str());

        if (status!=0){
            std::cout << "Something went wrong when attempting to use WFBP. Call exited with code: " << status << std::endl;
            std::cout << "Exiting." << std::endl;
            exit(1);
        }

        // Read the reconstruction from disk in allocated array
        std::string wfbp_recon_path=rp->output_dir + "/init_wfbp.img";

        std::cout << "WFBP PATH: " << wfbp_recon_path << std::endl;
        
        std::ifstream infile(wfbp_recon_path,std::ios::binary);
        if (!infile){
            std::cout << "File was not read properly." << std::endl;
            std::cout << "Could not initialize WFBP reconstruction volume from given file. Exiting" << std::endl;
            exit(1);
        }

        infile.read((char*)data->recon_volume,rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z*sizeof(float));
        infile.close();
        
        // Convert attenuation values from mm to cm
        for (int i=0; i<rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z; i++){
            data->recon_volume[i]=data->recon_volume[i]*10.0f;
        }

        // Flip each WFBP image if table direction is out (rp->table_direction>0)
        if (rp->table_direction>0){
            for (int k=0; k<rp->num_voxels_z; k++){
                for (int i=0; i<rp->num_voxels_x; i++){
                    for (int j=0; j<rp->num_voxels_y/2; j++){
                        size_t input_idx  = rp->num_voxels_x*j+i                      + rp->num_voxels_x*rp->num_voxels_y*k;
                        size_t output_idx = rp->num_voxels_x*(rp->num_voxels_y-1-j)+i + rp->num_voxels_x*rp->num_voxels_y*k;
        
                        float tmp=data->recon_volume[input_idx];
                        data->recon_volume[input_idx]=data->recon_volume[output_idx];
                        data->recon_volume[output_idx]=tmp;                
                    }            
                }
            }
        }

        // Transpose all WFBP slices to match coordinate frame of ICD
        float * tmp= new float[rp->num_voxels_x*rp->num_voxels_y];
        for (int k=0; k<rp->num_voxels_z; k++){
            for (int i=0; i<rp->num_voxels_x; i++){
                for (int j=0; j<rp->num_voxels_y; j++){                    
                    size_t input_idx   = rp->num_voxels_x*i+j + rp->num_voxels_x*rp->num_voxels_y*k;
                    size_t output_idx  = rp->num_voxels_x*j+i;
        
                    tmp[output_idx] = data->recon_volume[input_idx];
                }            
            }
            std::memcpy((char*)&data->recon_volume[k*rp->num_voxels_y*rp->num_voxels_x],tmp,sizeof(float)*rp->num_voxels_x*rp->num_voxels_y);
        }
        delete[] tmp;

        // Rotate the WFBP reconstruction to match what we'll have from ICD w/ rotating slices
        rotate_slices_fixed2rotating(rp,data);
        
        // Clean up reconstruction (We may not actually want to do this?)

    }

    /* No initialization of the reconstruction volume requested. Start at all 0s */
    else{ 
        // Matrix volume initialized to all zeros in setup.cpp
    }    
}

void build_wfbp_parameter_file(const struct recon_params * rp){

    std::string wfbp_parameter_filepath = rp->output_dir + "/init_wfbp.prm";

    std::ofstream wfbp_parameter_fid;

    wfbp_parameter_fid.open(wfbp_parameter_filepath);

    // Write all our parameters out to a file
    wfbp_parameter_fid << "Nrows: "                   << rp->Nrows                         << std::endl;
    wfbp_parameter_fid << "CollSlicewidth: "          << rp->CollSlicewidth                << std::endl;
    wfbp_parameter_fid << "StartPos: "                << rp->z_start*10.0                  << std::endl;
    wfbp_parameter_fid << "EndPos:   "                << rp->z_end*10.0                    << std::endl;
    wfbp_parameter_fid << "PitchValue: "              << rp->table_feed_per_rotation*10.0  << std::endl;
    wfbp_parameter_fid << "AcqFOV: "                  << rp->acquisition_fov*10.0          << std::endl;
    wfbp_parameter_fid << "ReconFOV: "                << rp->recon_fov*10.0                << std::endl;
    wfbp_parameter_fid << "Readings: "                << rp->Readings                      << std::endl;
    wfbp_parameter_fid << "Xorigin: "                 << "0.0"                             << std::endl;
    wfbp_parameter_fid << "Yorigin: "                 << "0.0"                             << std::endl;
    wfbp_parameter_fid << "TableVerticalOffset: "     << "0.0"                             << std::endl;
    wfbp_parameter_fid << "Zffs: "                    << rp->Zffs                          << std::endl;
    wfbp_parameter_fid << "Phiffs: "                  << rp->Phiffs                        << std::endl;
    wfbp_parameter_fid << "Scanner: "                 << "definitionas.scanner"            << std::endl;
    wfbp_parameter_fid << "FileType: "                << rp->FileType                      << std::endl;
    wfbp_parameter_fid << "FileSubType: "             << rp->FileSubType                   << std::endl;
    wfbp_parameter_fid << "RawOffset: "               << rp->RawOffset                     << std::endl;
    wfbp_parameter_fid << "Nx: "                      << rp->nx                            << std::endl;
    wfbp_parameter_fid << "Ny: "                      << rp->ny                            << std::endl;
    wfbp_parameter_fid << "ImageOrientationPatient: " << "N/A"                             << std::endl;
    wfbp_parameter_fid << "RawDataDir: "              << split_path(rp->sinogram_path)     << std::endl;
    wfbp_parameter_fid << "RawDataFile: "             << split_filename(rp->sinogram_path) << std::endl;
    wfbp_parameter_fid << "OutputDir: "               << rp->output_dir                    << std::endl;
    wfbp_parameter_fid << "OutputFile: "              << "init_wfbp.img"                   << std::endl;
    wfbp_parameter_fid << "ReconKernel: "             << "1"                               << std::endl;
    wfbp_parameter_fid << "SliceThickness: "          << rp->slice_thickness*10.0          << std::endl;
    wfbp_parameter_fid << "AdaptiveFiltration: "      << "1.0"                             << std::endl;
    /// done

    wfbp_parameter_fid.close();
}
