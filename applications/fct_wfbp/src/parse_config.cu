/* FreeCT_ICD is CT reconstruction Software */
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
#include "recon_structs.h"
#include "parse_config.h"
#include <yaml-cpp/yaml.h>

// Macro that actually reads our YAML and handles keyword detection
#define parse_item(TAG_NAME,TYPE) try{\
    std::cout << #TAG_NAME": " ;                                      \
    for (int i=0; i< 35-sizeof(#TAG_NAME);i++){std::cout << " ";}     \
    rp->TAG_NAME=config[#TAG_NAME].as<TYPE>();                        \
    std::cout << "FOUND ("  <<  rp->TAG_NAME << ")"  <<std::endl;           \
    }catch(YAML::RepresentationException& e){std::cout << "NOT FOUND" << std::endl;}

namespace{
  void parse_string(std::string tagname, std::string& value, YAML::Node& config){
    try{
      std::cout << tagname + ": " ;                                      
      for (int i=0; i< 35-tagname.length();i++){
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

void parse_config(std::string config_file, struct recon_params * rp){
    // Load our YAML config file
    std::cout << "Config File: " << config_file << std::endl;    
    YAML::Node config = YAML::LoadFile(config_file);
    std::cout << std::endl;

    //parse_item(raw_data_dir,std::string);
    std::string tmp;
    parse_string("raw_data_dir",tmp,config);
    strcpy(rp->raw_data_dir,tmp.c_str());
    parse_string("output_dir",tmp,config);
    strcpy(rp->output_dir,tmp.c_str());
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
    parse_item(tube_start_angle,float);
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
