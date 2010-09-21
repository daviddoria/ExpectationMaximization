cmake_minimum_required(VERSION 2.6)

PROJECT(ExpectationMaximization)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(ExpectationMaximization Example.cxx vtkExpectationMaximization.cxx Model.cxx)
TARGET_LINK_LIBRARIES(ExpectationMaximization vtkInfovis vtkViews vtkCharts)

