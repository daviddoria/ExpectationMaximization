# Allow headers in tests to be included like
# #include "ExpectationMaximization.h" rather than needing
# #include "ExpectationMaximization/ExpectationMaximization.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

# Test 1D
add_executable(ExpectationMaximization1DTest ExpectationMaximization1DTest.cpp)
target_link_libraries(ExpectationMaximization1DTest libExpectationMaximization)
add_test(ExpectationMaximization1DTest ExpectationMaximization1DTest)

# Test ND
add_executable(ExpectationMaximizationNDTest ExpectationMaximizationNDTest.cpp)
target_link_libraries(ExpectationMaximizationNDTest libExpectationMaximization)
add_test(ExpectationMaximizationNDTest ExpectationMaximizationNDTest)
