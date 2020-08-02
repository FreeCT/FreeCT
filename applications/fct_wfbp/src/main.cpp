#include <iostream>
#include <string>

#include <fct/FreeCTRead.h>
#include <parse_config.h>
#include <util.h>

#include <boost/filesystem.hpp>

#include <ReconstructionMachine.h>

#include <fstream>

void usage(){
  std::cout << "Usage: fct_wfbp [options] reconstruction_configuration.yaml" << std::endl;
  std::cout << "           -v: verbose" << std::endl;
  std::cout << "     -d <int>: CUDA device to utilize" << std::endl;
}

int main(int argc, char ** argv){

  std::string recon_config_filepath = "";
  bool flag_verbose = false;
  bool flag_testing = false;
  int cuda_device = 0;
  
  // Parse our command line inputs
  if (argc<2){
    usage();
  }
  
  for (int i=1; i<(argc-1);i++){
    if (i>=argc){
      break;
    }
    
    std::string arg = argv[i];
    
    if (arg=="-v")
      flag_verbose = true;
    else if (arg=="-t")
      flag_testing = true;
    else if (arg=="-d")
      cuda_device = std::stoi(argv[++i]);
    else{
      std::cout << "ERROR: Unrecognized argument \"" << arg << "\"" << std::endl;
      usage();
    }
  }

  recon_config_filepath = argv[argc-1];
  
  // Configure CUDA device
  bool was_successful = validate_selected_device(cuda_device);
  if (!was_successful){
    exit(1);
  }

  // Instantiate and set up our reconstruction
  fct::ReconstructionMachine rm;
  rm.LoadReconConfiguration(recon_config_filepath);
  rm.LoadRawData();
  rm.ConfigureCTGeometry();
  rm.PrintCTGeometry();
  rm.RunReconstruction();
  
  //////// Get our reconstruction parameters from the configuration file
  ////ReconConfig rp;
  ////if (!boost::filesystem::exists(recon_config_filepath)){
  ////  std::cout << "ERROR: Recon configuration filepath does not exist! (" << recon_config_filepath << ")" << std::endl;
  ////  exit(1);
  ////}
  ////
  ////parse_config(recon_config_filepath,rp);
  ////
  ////// Load our dataset
  ////std::string raw_data_path = rp.raw_data_dir;
  ////std::shared_ptr<fct::RawDataSet> ds = std::make_shared<fct::DicomDataSet>();
  ////ds->setPath(raw_data_path);
  ////ds->initialize();
  ////ds->readAll();
  ////
  ////// Create the CT geometry structure
  ////CTGeometry cg;
  ////configure_ct_geometry(ds,cg);
  ////
  ////// Get the data into CUDA-suitable shape (i.e. raw memory)
  ////std::shared_ptr<float> DATAPOINTER_raw(new float[cg.num_detector_rows*cg.num_detector_cols*cg.total_number_of_projections]);
  ////
  ////float * tmp_raw = DATAPOINTER_raw.get();
  ////for (size_t projection_idx=0; projection_idx < cg.total_number_of_projections;projection_idx++){
  ////  size_t data_offset = cg.num_detector_cols * cg.num_detector_rows * projection_idx;
  ////  ds->copyProjection(projection_idx,&tmp_raw[data_offset]);
  ////}
  ////
  ////if (flag_testing){
  ////  std::ofstream ofs("/home/john/Desktop/raw_debug.dat",std::ios::binary);
  ////  ofs.write((char*)tmp_raw,cg.num_detector_cols * cg.num_detector_rows * cg.total_number_of_projections*sizeof(float));
  ////}
  ////
  ////// Rebin and filter the data
  ////std::shared_ptr<float> DATAPOINTER_rebin(new float[cg.num_detector_rows*cg.num_detector_cols*cg.total_number_of_projections]);
  ////rebin(DATAPOINTER_rebin,DATAPOINTER_raw,cg,rp);
  ////
  ////float * tmp_rebin = DATAPOINTER_rebin.get();
  ////if (flag_testing){
  ////  std::ofstream ofs("/home/john/Desktop/rebin_debug.dat",std::ios::binary);
  ////  ofs.write((char*)tmp_rebin,cg.num_detector_cols * cg.num_detector_rows * cg.total_number_of_projections*sizeof(float));
  ////}
  
  return 0;
}