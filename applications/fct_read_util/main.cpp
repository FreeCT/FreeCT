#include <iostream>
#include <memory>

#include <FreeCTRead.h>

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
  if (output_dirpath==""){
    std::string message = "Output path must be set!";
    error_out(message);
  }

  // Ensure input and output directories/paths exist
  if (!boost::filesystem::exists(input_dirpath)){
    std::cout << "Could not find input dirpath: " << input_dirpath << std::endl;
    exit(1);
  }

  if (!boost::filesystem::exists(output_dirpath)){
    std::cout << "Could not find output dirpath: " << output_dirpath << std::endl;
    exit(1);
  }
   
  std::cout << "Input file/directory: " << input_dirpath << std::endl;
  std::cout << "Output directory:     " << output_dirpath << std::endl;
  
  // Runtime polymorphism to eventually support multiple raw data formats
  std::unique_ptr<fct::RawDataSet> ds = std::make_unique<fct::DicomDataSet>();

  //fct::RawDataSet * ds;
  //fct::DicomDataSet dicom_ds;
  //ds = &dicom_ds;

  ds->setPath(input_dirpath);
  ds->readAll();

  ds->writeAll(output_dirpath);

  return 0;
}