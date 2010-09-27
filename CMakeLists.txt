cmake_minimum_required(VERSION 2.6)

PROJECT(ExpectationMaximization)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(VXL REQUIRED)
INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)

ADD_EXECUTABLE(ExpectationMaximization Example.cxx vtkExpectationMaximization.cxx Gaussian1D.cxx GaussianND.cxx Model.cxx Plotting.cxx)
TARGET_LINK_LIBRARIES(ExpectationMaximization vtkInfovis vtkViews vtkCharts
vnl rrel #vgl vcl
)

