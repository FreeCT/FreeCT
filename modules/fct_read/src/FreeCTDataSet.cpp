#include "FreeCTDataSet.h"

#include <boost/filesystem.hpp>
#include <yaml-cpp/yaml.h>

namespace fct{

  bool FreeCTFrame::readFromFile(std::string filepath, size_t frame_idx){
    std::ifstream ifs(filepath,std::ios::binary);

    size_t offset = (3 + m_detector_rows*m_detector_channels)*sizeof(float)*frame_idx;

    ifs.seekg(offset,ifs.beg);
    ifs.read((char*)&m_dfc_angular_position,sizeof(float));
    ifs.read((char*)&m_dfc_axial_position,sizeof(float));
    ifs.read((char*)&m_dfc_radial_distance,sizeof(float));
    
    m_projection.resize(m_detector_rows,m_detector_channels);
    ifs.read((char*)m_projection.data(),m_detector_rows*m_detector_channels*sizeof(float));
    
    return true;
  }

  void FreeCTDataSet::readAll(){

    std::string meta_filepath = m_path + "/" + "meta.yaml";
    std::string source_filepath = m_path + "/" + "source_positions.dat";
    std::string projections_filepath = m_path + "/" + "projections.dat";

    // Ensure the necessary files exist
    if (!boost::filesystem::exists(meta_filepath)){
      std::cout << "ERROR: Could not find metadata filepath " + meta_filepath << std::endl;
      exit(1);
    }
    if (!boost::filesystem::exists(source_filepath)){
      std::cout << "ERROR: Could not find source positions filepath " + source_filepath << std::endl;
      exit(1);
    }
    if (!boost::filesystem::exists(projections_filepath)){
      std::cout << "ERROR: Could not find projection data filepath " + projections_filepath << std::endl;
      exit(1);
    }

    // Read and parse the meta.yaml file    
    YAML::Node doc = YAML::LoadFile(meta_filepath);
    doc >> *this;

    //m_data.resize(m_total_num_projections);
    
    // Read and parse the source positions and frame data
    for (size_t i=0;i<m_total_num_projections;i++){

      std::unique_ptr<fct::RawDataFrame> rdf = std::make_unique<fct::FreeCTFrame>();

      rdf->setDetectorRows(m_detector_rows);
      rdf->setDetectorChannels(m_detector_channels);
      rdf->readFromFile(m_path);

      m_data.push_back(std::move(rdf));
      
    }
    

  }
}