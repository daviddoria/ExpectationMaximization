cmake_minimum_required(VERSION 2.6)

PROJECT(ExpectationMaximization)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

if(NOT vtkViewsContext2D_LOADED)
  message(FATAL_ERROR "vtkViewsContext2D is required but not available")
endif()

if(NOT vtkChartsCore_LOADED)
  message(FATAL_ERROR "vtkChartsCore is required but not available")
endif()

FIND_PACKAGE(VXL REQUIRED)
INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)

add_subdirectory(KMeansClustering)

add_library(libExpectationMaximization vtkExpectationMaximization.cxx Gaussian1D.cxx GaussianND.cxx Model.cxx Plotting.cxx)
TARGET_LINK_LIBRARIES(libExpectationMaximization
libKMeansClustering
${VTK_LIBRARIES}
vnl vnl_algo)

ADD_EXECUTABLE(ExpectationMaximization Example.cxx)
TARGET_LINK_LIBRARIES(ExpectationMaximization
libExpectationMaximization
libKMeansClustering
${VTK_LIBRARIES}
vnl vnl_algo)

