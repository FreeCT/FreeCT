#include <iostream>
#include <memory>

#include <FreeCTRead.h>

int main(int argc, char ** argv){

  std::string test_path = std::string(argv[1]) + "/" + "freect_test";

  std::unique_ptr<fct::RawDataSet> rds = std::make_unique<fct::FreeCTDataSet>();
  rds->setPath(test_path);
  rds->readAll();
  rds->printMetadata();
  
}