add_executable(fct_read_util main.cpp)
target_include_directories(fct_read PUBLIC ${INCLUDE_DIRS})
target_link_libraries(fct_read_util fct_read Boost::system Boost::filesystem ${DCMTK_LIBRARIES} yaml-cpp)
install(TARGETS fct_read_util DESTINATION bin/)