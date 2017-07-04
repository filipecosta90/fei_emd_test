cmake_minimum_required (VERSION 2.8)
project (hdf_wrapper_test)

include_directories("${PROJECT_BINARY_DIR}/../..")

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
link_directories( ${HDF5_LIBRARY_DIRS} )
include_directories( ${HDF5_INCLUDE_DIRS} )

add_executable(emd_test test_emd.cpp)
target_link_libraries(emd_test ${HDF5_CXX_LIBRARIES} ) 
