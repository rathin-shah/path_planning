cmake_minimum_required(VERSION 3.0.2)
project(interview-planning)

set(CMAKE_CXX_STANDARD 17)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

find_package(Eigen3 REQUIRED)

add_executable(CheckPath src/main.cpp)
target_include_directories(CheckPath PRIVATE include)
target_link_libraries(CheckPath PRIVATE Eigen3::Eigen)
