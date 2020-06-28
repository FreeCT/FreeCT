cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
set(PROJECT_NAME FreeCT LANGUAGES CXX CUDA)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-pthread -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -g")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Check dependencies
FIND_PACKAGE(Boost COMPONENTS system filesystem REQUIRED)
IF (Boost_FOUND)
  message("FOUND BOOST!")
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
  ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

find_package(yaml-cpp REQUIRED)
if (yaml-cpp_FOUND)
  message("Found yaml-cpp.  Assuming it's installed on library search path.")
ELSE(yaml-cpp_FOUND)
endif(yaml-cpp_FOUND)

FIND_PACKAGE(DCMTK)
if (DCMTK_FOUND)
  message("FOUND DCMTK!")
ELSE(DCMTK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot find DCMTK, are you sure you've installed it?")
ENDIF(DCMTK_FOUND)

set(INCLUDE_DIRS modules/fct_read/include/ ${DCMTK_INCLUDE_DIRS})

add_subdirectory(modules)
add_subdirectory(applications)
