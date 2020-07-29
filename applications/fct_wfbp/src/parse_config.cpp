#include <parse_config.h>
#include <iostream>
#include <yaml-cpp/yaml.h>

// Macro that actually reads our YAML and handles keyword detection
// I don't love doing this with a macro, but it works pretty well
// so I'm not going to change it for the time being.
#define parse_item(TAG_NAME,TYPE) try{                                  \
    std::cout << #TAG_NAME": " ;                                        \
    for (size_t i=0; i< (size_t)35-sizeof(#TAG_NAME);i++){std::cout << " ";}  \
    rp.TAG_NAME=config[#TAG_NAME].as<TYPE>();                          \
    std::cout << "FOUND ("  <<  rp.TAG_NAME << ")"  <<std::endl;       \
  }catch(YAML::RepresentationException& e){std::cout << "NOT FOUND" << std::endl;}

namespace{
  void parse_string(std::string tagname, std::string& value, YAML::Node& config){
    try{
      std::cout << tagname + ": " ;                                      
      for (size_t i=0; i< 35-tagname.length()-1;i++){
        std::cout << " ";
      }     
      value=config[tagname].as<std::string>();
      std::cout << "FOUND ("  <<  value << ")"  << std::endl;
    }
    catch(YAML::RepresentationException& e){
      std::cout << "NOT FOUND" << std::endl;
    }    
  }
}

void parse_config(std::string config_file, ReconConfig& rp){
  // Load our YAML config file
  std::cout << "Config File: " << config_file << std::endl;    
  YAML::Node config = YAML::LoadFile(config_file);
  std::cout << std::endl;

  //parse_item(raw_data_dir,std::string);rp- std::string tmp;
  std::string tmp;
  parse_string("raw_data_dir",tmp,config);
  strcpy(rp.raw_data_dir,tmp.c_str());
  parse_string("output_dir",tmp,config);
  strcpy(rp.output_dir,tmp.c_str());
  strcpy(rp.output_file,(config_file+"_recon.dat").c_str());
  //parse_item(output_dir,std::string);
    
  parse_item(start_pos,float);
  parse_item(end_pos,float);
  parse_item(recon_fov,double);
  parse_item(slice_thickness,double);
  parse_item(nx,size_t);
  parse_item(ny,size_t);
  parse_item(recon_kernel,int);
  parse_item(x_origin,float);
  parse_item(y_origin,float);
  parse_item(tube_angle_offset,float);
  parse_item(adaptive_filtration_s,float);
    
  // Deprecated
  // ========================================
  // parse_item(dx,double);
  // parse_item(dy,double);
  // parse_item(dz,double);
  // parse_item(fov_radius,double);
  // parse_item(center_voxel_x,double);
  // parse_item(center_voxel_y,double);
  // parse_item(center_voxel_z,double);
  // parse_item(nz,size_t);
}

void configure_ct_geometry(std::shared_ptr<fct::RawDataSet> ds,CTGeometry& cg){
  // Physical geometry of the scanner (cannot change from scan to scan)
  
  cg.total_number_of_projections  = ds->getTotalNumProjections();
  cg.projections_per_rotation     = ds->getProjectionsPerRotation();
  cg.detector_pixel_size_col     = ds->getDetectorTransverseSpacing();
  cg.detector_pixel_size_row     = ds->getDetectorAxialSpacing();
  cg.num_detector_cols            = ds->getDetectorChannels();
  cg.num_detector_rows            = ds->getDetectorRows();
  cg.detector_central_col         = ds->getDetectorCentralChannel();
  cg.detector_central_row         = ds->getDetectorCentralRow();
  cg.distance_source_to_isocenter = ds->getDistSourceToIsocenter(); 
  cg.distance_source_to_detector  = ds->getDistSourceToDetector();
    
  cg.collimated_slice_width       = ds->getDistSourceToIsocenter()*(ds->getDetectorAxialSpacing()/ds->getDistSourceToDetector());
  
  //cg.z_rot = ds->getTablePosition(cg.projections_per_rotation-1) - ds->getTablePosition(0);
  cg.z_rot = fabs(ds->getTablePosition(cg.projections_per_rotation) - ds->getTablePosition(0));
  
  float detector_cone_offset = ((float)(cg.num_detector_rows - 1))/2.0f; // May not be 100% accurate if central detector is not necessarily in the middle
  cg.theta_cone=2.0f*atan(detector_cone_offset * cg.collimated_slice_width/cg.distance_source_to_isocenter);

  cg.acquisition_field_of_view = 2.0f * cg.distance_source_to_isocenter*sin((float(cg.num_detector_cols-1.0f)/2.0f) * ds->getDetectorTransverseSpacing() * (1.0f/cg.distance_source_to_detector));
  
  std::cout << "CT Geometry and Scan derived parameters: "   << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "Num projections per turn:           "        << cg.projections_per_rotation     << std::endl;
  std::cout << "Num detector channels:              "        << cg.num_detector_cols            << std::endl;
  std::cout << "Num detector rows:                  "        << cg.num_detector_rows            << std::endl;
  std::cout << "Radius src->isocenter (mm):         "        << cg.distance_source_to_isocenter << std::endl;
  std::cout << "Radius src->detector (mm):          "        << cg.distance_source_to_detector  << std::endl;
  std::cout << "Table feed per rotation (mm):       "        << cg.z_rot                        << std::endl;
  std::cout << "Theta cone (rad):                   "        << cg.theta_cone                   << std::endl;
  std::cout << "Central channel:                    "        << cg.detector_central_col         << std::endl;
  std::cout << "Central row:                        "        << cg.detector_central_row         << std::endl;
  std::cout << "Acquisition FOV (mm):               "        << cg.acquisition_field_of_view    << std::endl;  
  std::cout << "Collimated slicewidth at isocenter: "        << cg.collimated_slice_width       << std::endl;
  
}