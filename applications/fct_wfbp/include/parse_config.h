#pragma once

#include <recon_structs.h>
#include <string>
#include <cstring>
#include <memory>
#include <fct/FreeCTRead.h>

void parse_config(std::string configuration_file, ReconConfig& rp);
void configure_ct_geometry(std::shared_ptr<fct::RawDataSet> ds,CTGeometry& cg);