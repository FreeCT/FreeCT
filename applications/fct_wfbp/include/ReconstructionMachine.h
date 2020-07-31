#include <fct/FreeCTRead.h>
#include <recon_structs.h>
#include <memory>
#include <cstring>

namespace fct{
  
  class ReconstructionMachine{
  public:
    ReconstructionMachine(){
      memset(&m_rp,0,sizeof(ReconConfig));
      memset(&m_cg,0,sizeof(CTGeometry));
    };
    ~ReconstructionMachine(){};

    void LoadReconConfiguration(std::string filepath);
    void LoadRawData();
    void ConfigureCTGeometry();
    void RunReconstruction();
    void PrintCTGeometry();

  private:

    void AllocateKeyArrays();
    void RebinAndFilter();
    void Backproject();
    void ApplyFinalSliceThicknessing();

    ReconConfig m_rp;
    CTGeometry m_cg;

    // Key host pointers (data is on HOST)
    std::shared_ptr<fct::RawDataSet> m_org_data_set;
    std::shared_ptr<float> m_raw_data;
    std::shared_ptr<float> m_reconstructed_data;
    
    std::vector<float> m_tube_angles;
    std::vector<float> m_table_positions;
    std::vector<float> m_slice_locations_collimated_slice_width;
    std::vector<float> m_slice_locations_requested_slice_width;

    // Key device pointers (data is on GPU)
    // Data here is shared between multiple methods facilitating the reconstruction
    float * m_d_filtered_projection_data;
    float * m_d_reconstruction_collimated_slice_width;
    float * m_d_slice_locations_collimated_slice_width;
  };

}