cmake_minimum_required (VERSION 2.8)
project (hdf_wrapper_test)

include_directories("${PROJECT_BINARY_DIR}/../..")

find_package(OpenCV  3.0.0 REQUIRED)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
find_package(Boost 1.64.0 COMPONENTS log log_setup system program_options filesystem thread iostreams REQUIRED)

link_directories( ${HDF5_LIBRARY_DIRS} )
include_directories( ${HDF5_INCLUDE_DIRS} )
include_directories( ${OpenCV_INCLUDE_DIRS} )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories( ${Boost_LIBRARY_PATH} )

add_executable(emd_test test_emd.cpp)
target_link_libraries(emd_test ${HDF5_CXX_LIBRARIES}  ${Boost_LIBRARIES} ${OpenCV_LIBS} )
