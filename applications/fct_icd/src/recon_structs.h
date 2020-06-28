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

#pragma once

#include <string>

struct recon_params {

    // Pertinent paths
    std::string initial_recon_path;     // Full path to initializer recon
    std::string matrix_path;            // Full path to system matrix
    std::string output_dir;             // Full path to output directory
    std::string output_file;            // Final reconstruction filename (e.g. my_reconstruction.img)
    std::string potential;              // Potential function to use for penalty term
    std::string recon_path;             // full path to output file? (I have no idea what the difference is here from output_dir)
    std::string sinogram_path;          // full path to raw data file

    // Scanner Geometry
    double acquisition_fov;             // cm (scanner maximum FOV)
    double anode_angle;                 // degrees
    double axial_detector_spacing;      // cm 
    double axial_focal_spot_shift;      // cm 
    double center_channel_non_ffs;      // voxels
    double center_row;                  // voxels
    double focal_spot_radius;           // cm 
    size_t n_channels; 
    size_t num_views_per_turn;
    size_t num_views_per_turn_without_ffs;
    double source_detector_distance;    // cm    
    double transaxial_detector_spacing; // cm 
    double transaxial_focal_spot_shift; // cm
    double fan_angle_increment;         // radians

    // Iterative recon parameters
    std::string penalty;                // Either "quadratic" or "edge-preserving"
    double lambda;                      // Regularizer parameter
    double delta;                       // Edge-preserving parameter. (0.005 is a good value).
    size_t num_iterations;

    // Scan specific parameters
    size_t Nrows;                       //
    double first_view_angle;            // degrees
    size_t num_z_ffs;
    size_t num_phi_ffs;
    size_t num_views;    
    double table_feed_per_rotation;     // cm
    double tube_angle_increment;        // radians
    double z_end;                       //cm 
    double z_first_view;                // cm 
    double z_good_slice;                // cm 
    double z_matrix_start;              // cm 
    double z_start;                     // cm
    int table_direction;                // -1 if into the scanner +1 if out of the scanner
    
    // Recon Geometry
    double fov_radius;                  // Reconstruction FOV (default value 50.0 cm)
    double center_voxel_x;              // voxels
    double center_voxel_y;              // voxels 
    double center_voxel_z;              // voxels 
    double dx;                          // cm voxel size
    double dy;                          // cm voxel size
    double dz;                          // cm voxel size
    size_t nx;                          // voxels
    size_t ny;                          // voxels 
    size_t nz;                          // voxels

    // John's Additions
    size_t FileType;
    size_t FileSubType;
    size_t RawOffset;
    size_t Readings;
    size_t Zffs;
    size_t Phiffs;
    double CollSlicewidth;
    size_t wfbp_initialize;

    // Derived values
    size_t Nrows_projection; // This is our adjusted value to account for z flying focal spot
    size_t views_per_slice;
    double voxel_size_x;
    double voxel_size_y;
    double voxel_size_z;
    double tube_z_increment;
    double beam_width;
    double beam_width_at_roi_edge;
    size_t num_views_for_system_matrix;
    
    size_t num_voxels_x;
    size_t num_voxels_y;
    size_t num_voxels_z;

    // More updates
    double slice_thickness;
    double recon_fov;
    
};

struct ct_data{
    float * system_matrix;
    float * raw;
    float * recon_volume;
    float * final_image_stack;
    double * tube_angles;
    double * table_positions;
    double * slice_locations;
    size_t * slice_indices;
};
