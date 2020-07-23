#include <iostream>
#include <memory>

#include <fct/FreeCTRead.h>

#include <boost/filesystem.hpp>

void usage(){
  printf("freect_read -i /path/to/dicom_directory/ -o /path/to/output_directory/\n");
  printf("Copyright (C) John Hoffman 2020\n");
};

void error_out(std::string message){
  std::cout << message << std::endl;
  usage();
  exit(1);
}

int main(int argc, char ** argv){
  
  std::string input_dirpath = "";
  std::string output_dirpath = "";
  std::string recon_filepath = "";
  bool single_slice_mode = false;
  
  // Parse command line arguments
  if (argc<5){
    std::cout << "Not enough input arguments!" << std::endl;
    usage();
    exit(1);
  }
  
  for (int i=1;i<argc;i++){

    std::string arg(argv[i]);
      
    if (arg=="-i")
      input_dirpath = argv[++i];
      
    else if (arg=="-o")
      output_dirpath = argv[++i];

    else if (arg=="-r")
      recon_filepath = argv[++i];

    else if (arg=="--single-slice")
      single_slice_mode = true;
    
    else{
      std::string message = "Unrecognized option \"" + std::string(argv[i]) +  "\" requested";
      error_out(message);
    }
  }

  // Validate the required inputs were set
  if (input_dirpath==""){
    std::string message = "Input path must be set!";
    error_out(message);
  }
  if (output_dirpath=="" && recon_filepath==""){
    std::string message = "Output path or recon_filepath must be set!";
    error_out(message);
  }

  // Ensure input and output directories/paths exist
  if (!boost::filesystem::exists(input_dirpath)){
    std::cout << "Could not find input dirpath: " << input_dirpath << std::endl;
    exit(1);
  }

  if (output_dirpath!="" && !boost::filesystem::exists(output_dirpath)){
    std::cout << "Could not find output dirpath: " << output_dirpath << std::endl;
    exit(1);
  }
   
  std::cout << "Input file/directory: " << input_dirpath << std::endl;
  std::cout << "Output directory:     " << output_dirpath << std::endl;
  
  // Runtime polymorphism to eventually support multiple raw data formats
  std::unique_ptr<fct::RawDataSet> ds = std::make_unique<fct::DicomDataSet>();
  ds->setPath(input_dirpath);
  ds->initialize();
  ds->readAll();

  if (recon_filepath!="")
    ds->writeReconFile(recon_filepath,single_slice_mode);
  
  if (output_dirpath!="")
    ds->writeAll(output_dirpath);

  return 0;
}