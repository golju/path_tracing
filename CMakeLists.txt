cmake_minimum_required(VERSION 3.17)
project(path_tracing)

set(CMAKE_CXX_STANDARD 17)
set(LIBRARY_NAME librender)

set(${PROJECT_NAME}_SOURCE_FILES
        RenderLibrary.h
        RenderLibrary.cpp
        )

# Packages
find_package(OpenCV REQUIRED)
include_directories( ${OpenCV_INCLUDE_DIRS} )

# New library
add_library(${LIBRARY_NAME} SHARED ${${PROJECT_NAME}_SOURCE_FILES} )

add_executable(path_tracing main.cpp ${${PROJECT_NAME}_SOURCE_FILES})
target_link_libraries(path_tracing PRIVATE opencv_core opencv_highgui ${LIBRARY_NAME})