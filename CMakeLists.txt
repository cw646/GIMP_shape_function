project(gimp_test)

set(CMAKE_CXX_STANDARD 14)

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(gimp_test main.cpp)
