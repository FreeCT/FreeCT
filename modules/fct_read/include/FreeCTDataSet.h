#include "RawDataSet.h"

namespace fct{
  
  class FreeCTFrame: public RawDataFrame{
    friend class FreeCTDataSet;
    virtual bool readFromFile(std::string filepath, size_t frame_idx);
  };
  
  class FreeCTDataSet: public RawDataSet{
    
    virtual void readAll();
  };
  
}
