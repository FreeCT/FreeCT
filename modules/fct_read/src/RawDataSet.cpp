#include "RawDataSet.h"

#include <yaml-cpp/yaml.h>

#include <cstring>
#include <fstream>

#include <boost/filesystem.hpp>

namespace fct{

  /*  METHODS FOR RawDataSet */
  void RawDataSet::setPath(std::string path){
    m_path = path;
  }

  void RawDataSet::getYAMLEmitter(YAML::Emitter& ye){
    ye << YAML::BeginMap;

    ye << YAML::Comment("FREECT READER OUTPUT FOR: " + m_path) << YAML::Newline << YAML::Newline;
    ye << YAML::Comment("DETECTOR");
    ye << YAML::Key << "manufacturer" << YAML::Value << m_manufacturer;
    ye << YAML::Key << "detector_rows" << YAML::Value << m_detector_rows;
    ye << YAML::Key << "detector_channels" << YAML::Value << m_detector_channels;
    ye << YAML::Key << "detector_transverse_spacing" << YAML::Value << m_detector_transverse_spacing;
    ye << YAML::Key << "detector_axial_spacing" << YAML::Value << m_detector_axial_spacing;
    ye << YAML::Key << "detector_shape" << YAML::Value << m_detector_shape;
    
    ye << YAML::Newline << YAML::Newline;
    ye << YAML::Comment("SCANNER GEOMETRY");
    ye << YAML::Key << "distance_source_to_detector" << YAML::Value << m_dist_source_to_detector;
    ye << YAML::Key << "distance_source_to_isocenter" << YAML::Value << m_dist_source_to_isocenter;
    ye << YAML::Key << "detector_central_row" << m_detector_central_row;
    ye << YAML::Key << "detector_central_channel" << m_detector_central_channel;
    ye << YAML::Key << "projection_geometry" << YAML::Value << m_projection_geometry;
    
    ye << YAML::Newline << YAML::Newline;
    ye << YAML::Comment("SCAN");
    ye << YAML::Key << "scan_type" <<  YAML::Value << m_scan_type;
    ye << YAML::Key << "projections_per_rotation" <<  YAML::Value << m_projections_per_rotation;
    ye << YAML::Key << "flying_focal_spot_mode" <<  YAML::Value << m_flying_focal_spot_mode;
    
    ye << YAML::Key << "total_num_projections" <<  YAML::Value << m_total_num_projections;
    ye << YAML::EndMap;    
  };

  void RawDataSet::writeAll(std::string dirpath){
    // Emit YAML metadata and save to disk
    YAML::Emitter ye;
    getYAMLEmitter(ye);

    std::ofstream ofs_meta(dirpath + "/" + "meta.yaml");
    ofs_meta << ye.c_str() << std::endl;

    std::cout << ye.c_str() << std::endl;

    // Save source positions and projection data
    std::ofstream out(dirpath + "/" + "projections.dat", std::ios::out | std::ios::binary | std::ios::trunc);
    
    for (auto &f: m_data){
      float  angle = f->getDFCAngularPosition() + f->getFFSAngularShift();
      float  axial = f->getDFCAxialPosition() + f->getFFSAxialShift();
      float radial = f->getDFCRadialPosition() + f->getFFSRadialShift();
      out.write((char*)&angle ,sizeof(float));
      out.write((char*)&axial ,sizeof(float));
      out.write((char*)&radial,sizeof(float));
      out.write((char*)f->m_projection.data(),m_detector_rows*m_detector_channels*sizeof(float));
    }
  }

  void RawDataSet::printMetadata(){
    YAML::Emitter ye;
    getYAMLEmitter(ye);
    std::cout << ye.c_str() << std::endl;
  }


  /* METHODS FOR RawDataFrame */
  int RawDataFrame::getRows(){
    return m_detector_rows;
  };
  
  int RawDataFrame::getCols(){
    return m_detector_channels;
  };

  
}