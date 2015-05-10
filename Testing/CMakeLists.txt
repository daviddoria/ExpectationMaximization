# Allow headers in tests to be included like
# #include "ExpectationMaximization.h" rather than needing
# #include "ExpectationMaximization/ExpectationMaximization.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

# Test 1D
add_executable(ExpectationMaximization1DTest ExpectationMaximization1DTest.cpp)
target_link_libraries(ExpectationMaximization1DTest libExpectationMaximization)
add_test(ExpectationMaximization1DTest ExpectationMaximization1DTest)

# Test ND
add_executable(ExpectationMaximization2DTest ExpectationMaximization2DTest.cpp)
target_link_libraries(ExpectationMaximization2DTest libExpectationMaximization)
add_test(ExpectationMaximization2DTest ExpectationMaximization2DTest)
