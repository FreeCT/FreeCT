add_executable(FreeCTDataSet_test FreeCTDataSet_test.cpp)
target_include_directories(FreeCTDataSet_test PUBLIC ${INCLUDE_DIRS})
target_link_libraries(FreeCTDataSet_test fct_read Boost::system Boost::filesystem ${DCMTK_LIBRARIES} yaml-cpp)