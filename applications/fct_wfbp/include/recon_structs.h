#pragma once

#include <cstddef>
#include <cmath>

struct CTGeometry{
  size_t total_number_of_projections; //
  size_t projections_per_rotation; //
  float detector_pixel_size_col;//
  float detector_pixel_size_row;//
  size_t num_detector_cols; // "Channels" //
  size_t num_detector_rows;//
  float detector_central_col;//
  float detector_central_row;//
  float distance_source_to_detector;//
  float distance_source_to_isocenter;//
  float collimated_slice_width; // Native detector slice thickness
  float theta_cone;//
  float acquisition_field_of_view;//
  float z_rot; // table feed per rotation

  // Not handled yet, but may need to be in the future
  float anode_angle;//cg.anode_angle=7.0f*pi/180.0f;
  bool ffs_phi;
  bool ffs_z;
  bool ffs_diag;

};

struct ReconConfig{
  char raw_data_dir[4096];
  char output_dir[4096];
  char output_file[255];

  float start_pos;
  float end_pos;
  float recon_fov;
  float slice_thickness;
  float slice_pitch;
  size_t nx;
  size_t ny;
  int recon_kernel;
  float x_origin;
  float y_origin;
  float tube_angle_offset;
  float adaptive_filtration_s;
};

struct GPUPrecompute{

  float half_acquisition_fov;
  float half_acquisition_fov_squared;
  int projections_per_half_turn;
  int n_half_turns;
  float distance_source_to_isocenter_squared;
  float z_rot_over_2_pi;
  float recip_distance_source_to_detector;
  float recip_distance_source_to_isocenter;
  float tanf_theta_cone;
  float pixel_scale;

  int n_slices_native;
  int n_slices_requested;

  void InitFromCTGeometry(CTGeometry cg){

    half_acquisition_fov                 = 0.5f*cg.acquisition_field_of_view;
    half_acquisition_fov_squared         = half_acquisition_fov * half_acquisition_fov;
    projections_per_half_turn            = cg.projections_per_rotation/2;
    n_half_turns                         = floor(cg.total_number_of_projections/projections_per_half_turn);
    distance_source_to_isocenter_squared = cg.distance_source_to_isocenter*cg.distance_source_to_isocenter;
    z_rot_over_2_pi                      = cg.z_rot/(2.0f*3.1415926535897f);
    recip_distance_source_to_detector    = 1.0f/cg.distance_source_to_detector;
    recip_distance_source_to_isocenter   = 1.0f/cg.distance_source_to_isocenter;
    tanf_theta_cone                      = tan(cg.theta_cone/2.0f);
    pixel_scale                          = cg.distance_source_to_detector/(cg.distance_source_to_isocenter*cg.detector_pixel_size_col);
  }
  
};