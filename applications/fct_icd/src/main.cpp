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
#include <string>

#include "recon_structs.h"
#include "setup.h"
#include "initialize_recon_volume.h"
#include "generate_system_matrix.h"
#include "generate_system_matrix_ffs.h"
#include "rotate_slices.h"
#include "icd_iteration.h"

struct flags {
    bool testing;
    bool verbose;
    bool timing;
};

void empty_parameter_file();

void usage(){
    printf("\n");
    printf("usage: ctbb_icd [options] input_prm_file\n\n");
    printf("    Options:\n");
    printf("          -v: verbose.\n");
    printf(" --empty-prm: Print empty parameter file to stdout. (redirect output into .yaml file to store)\n");
    printf("    --timing: Display timing information for each step of the recon process\n");
    printf("\n");
    printf("Copyright Stefano Young, John Hoffman 2017\n\n");
    exit(0);
}

inline bool exists(const std::string& name){
    std::ifstream f(name.c_str());
    return f.good();
}

int main(int argc, char ** argv){

    struct ct_data data={};
    struct recon_params rp={}; // initialize parameter structure to zero
    struct flags flags={};
    
    /*--- Parse command line options ---*/ 
    if (argc<2)
        usage();

    if (argc==2){
        std::string curr_arg = argv[1];
        if (curr_arg.compare("--empty-prm")==0)
            empty_parameter_file();            
    }

    for (int i=1; i<argc-1; i++){
        std::string curr_arg = argv[i];

        if (curr_arg.compare("-t")==0)
            flags.testing=true;
        else if (curr_arg.compare("-v")==0)
            flags.verbose=true;
        else if (curr_arg.compare("--timing")==0)
            flags.timing=true;
        else
            usage();
    }
    
    /*--- Parse our configuration file (parse_config.cpp)---*/
    std::string parameter_file = argv[argc-1];
    rp=configure_recon_params(parameter_file);

    /*--- Get raw data (setup.cpp) ---*/
    load_raw(&rp,&data);

    /*--- All parameters now set to final values, map to CONSTANT ---*/
    const struct recon_params rp_const=rp;
    
    /*--- Initialize reconstruction volume (setup.cpp) ---*/
    // Perform wFBP reconstruction if using as input to ICD
    initialize_recon_volume(&rp_const,&data);

    /*--- Generate system matrix (generate_system_matrix.cpp or generate_system_matrix_ffs.cpp) ---*/
    // If matrix file does not exist, generate
    if (!exists(rp_const.matrix_path)){
        std::cout << "No existing matrix file found for reconstruction." << std::endl;
        if (!rp_const.Zffs && !rp_const.Phiffs)
            generate_system_matrix(&rp_const,&data);
        else
            generate_system_matrix_ffs(&rp_const,&data);
    }
    else{
        std::cout << "Existing matrix file FOUND." << std::endl;
        std::cout << "Using matrix file: " << rp_const.matrix_path << std::endl;
    }

    std::ofstream debug_outfile("/home/john/Desktop/debug_file.bin",std::ios_base::binary | std::ios::out);
    std::cout << "Number of elements: " << rp.nx*rp.ny*rp.num_voxels_z << std::endl;
    debug_outfile.write((char*)&data.recon_volume[0],rp.nx*rp.ny*rp.num_voxels_z*sizeof(float));
    debug_outfile.close();
    std::cout << "Debug file written to desktop." << std::endl;

    /*--- Perform ICD iterations (icd_iteration.cpp) ---*/
    icd_iteration(&rp_const,&data);

    /*--- De-rotate our ICD slices (rotate_slices.cpp) ---*/
    rotate_slices_rotating2fixed(&rp_const,&data);

    /*--- Write data to disk and clean up (setup.cpp) ---*/
    clean_up(&rp,&data);
    
    return 0;    
}

void empty_parameter_file(){
    std::cout << ""
        "# Empty parameter file for a FreeCT_ICD reconstruction\n"
        "# All units are CM, Radians, pixel units (unitless)\n"
        "# Config file syntax is YAML (http://yaml.org/) (i.e. no tabs)\n"
        "\n"
        "# Paths\n"
        "sinogram_path:      # Required, full path to raw projection data file\n"
        "output_dir:         # Required, full path to output directory\n"
        "output_file:        # Required, filename for final reconstruction file\n"
        "initial_recon_path: # Optional, full path to reconstruction file\n"
        "matrix_path:        # Optional, full path to existing stored matrix file\n"
        "\n"
        "# Raw data file information\n"
        "FileType:        # Required, choices are 0,1, or 4. 0: binary file; 1: Siemens PTR file; 4: Siemens IMA file\n"
        "FileSubType:     # Required, (if using IMA format). Specifies format of IMA-wrapped file.\n"
        "RawOffset:       # Required, offset in bytes to beginning of raw projection data (Note, this is often just 0)\n"
        "Readings:        # Required, total number of projections in raw data file to be reconstructed\n"
        "\n"
        "# Scanner Geometry\n"
        "# Parameters here are FIXED and should not change from one reconstruction to another\n"
        "acquisition_fov:                # Required, scanner FOV in cm\n"
        "n_channels:                     # Required, number of detector channels\n"
        "num_views_per_turn_without_ffs: # Required, number of projections per gantry rotation\n"
        "focal_spot_radius:              # Required, distance from focal spot to isocenter\n"
        "source_detector_distance:       # Required, distance from focal spot to detector center\n"
        "anode_angle:                    # Required (for flying focal spot scans), anode angle \n"
        "axial_detector_spacing:         # Required, detector spacing in CM AT THE DETECTOR in the Z direction (i.e. row separation)\n"
        "axial_focal_spot_shift:         # Required (for flying focal spot scans), Z direction shift in CM of focal spot\n"
        "center_channel_non_ffs:         # Required, central detector index (accounting for any quarter detector offset, if present)\n"
        "center_row:                     # Required, central row index (where line connecting source and isocenter intercepts detector plane)\n"
        "transaxial_detector_spacing:    # Required, detector spacing in CM AT THE DETECTOR in the axial plane (i.e. channel separation)\n"
        "transaxial_focal_spot_shift:    # Required (for flying focal spot scans), x-y direction shift in CM of focal spot\n"
        "\n"
        "# Scan specific parameters\n"
        "table_feed_per_rotation:    # Required, table feed in CM per rotation\n"
        "Zffs:                       # Required, 0 if no Z flying focal spot, 1 if Z flying focal spot used (Note: as of 2017-08-18, FreeCT_ICD does not reconstruct data with flying focal spots)\n"
        "Phiffs:                     # Required, 0 if no Phi (in-plane) flying focal spot, 1 if Phi flying focal spot used (Note: as of 2017-08-18, FreeCT_ICD does not reconstruct data with flying focal spots)\n"
        "Nrows:                      # Required, width of detector row at isocenter (note, this is the first number given in the collimation i.e. 16x1.2)\n"
        "CollSlicewidth:             # Required, width of detector row at isocenter (note, this is the second number given in the collimation i.e. 16x1.2)\n"
        "\n"
        "# Recon Geometry\n"
        "recon_fov:       # Required, field of view in CM of the reconstruction (i.e. reconstruction diameter)\n"
        "nx:              # Required, number of pixels in x-direction (i.e. width)\n"
        "ny:              # Required, number of pixels in y-direction (i.e. height)\n"
        "slice_thickness: # Required, slice thickness in CM (NOTE: Due to rotation slices approach, not all slice thickness are possible. See Xu et al. Phys. Med. Biol., 2012)\n"
        "z_start:         # Required, slice location of the first slice in the reconstruction stack\n"
        "z_end:           # Required, slice location of the last slice in the reconstruction stack\n"
        "\n"
        "# Iterative recon parameters\n"
        "wfbp_initialize:        # Required, 0 or 1. 1: Perform a WFBP reconstruction prior to running ICD iterations. 0: initialize reconstruction from all zeros.\n"
        "penalty:                # Required, string specifying penalty function to use. Choices are quadratic (default) or edge-preserving\n"
        "lambda:                 # Required, iteration reconstruction parameter. Default value is 0.1\n"
        "delta:                  # Required, edge-preserving parameter. Default value is 0.005\n"
        "system_memory_fraction: # Currently not used. May be utilized in a future release\n"
        "num_iterations:         # Required, number of ICD iteration to perform. W/ wFBP initialization, 50 is common. W/O wFBP initialization, 100-150 iterations typically required\n"
        "" << std::endl;
    exit(0);
}
