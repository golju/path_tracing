cmake_minimum_required(VERSION 3.17)
project(path_tracing)

set(CMAKE_CXX_STANDARD 17)

# Packages
find_package(OpenCV REQUIRED)
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(path_tracing main.cpp)

target_link_libraries( path_tracing PRIVATE opencv_core opencv_highgui )