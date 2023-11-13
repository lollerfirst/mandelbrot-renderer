cmake_minimum_required(VERSION 3.20)
project(mandelbrot-renderer CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP REQUIRED)
find_package(PNG REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(mandelbrot mandelbrot.cpp)

target_link_libraries(mandelbrot PRIVATE OpenMP::OpenMP_CXX PNG::PNG)
