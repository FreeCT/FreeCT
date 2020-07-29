#pragma once

#include <cstddef>

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
  size_t nx;
  size_t ny;
  int recon_kernel;
  float x_origin;
  float y_origin;
  float tube_angle_offset;
  float adaptive_filtration_s;
};