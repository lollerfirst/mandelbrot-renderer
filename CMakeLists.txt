cmake_minimum_required(VERSION 3.20)
project(mandelbrot-renderer CXX)

set(CMAKE_CXX_STANDARD 20)

if(APPLE)
    set(OPENMP_LIBRARIES "/usr/local/opt/libomp/lib")
    set(OPENMP_INCLUDES "/usr/local/opt/libomp/include")
    set(OpenMP_CXX_FLAGS "-L ${OPENMP_LIBRARIES} -lomp")
    set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS})
    include_directories(mandelbrot-renderer ${OPENMP_INCLUDES})
endif()

find_package(OpenMP REQUIRED)
find_package(PNG REQUIRED)
find_package(Boost REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(mandelbrot mandelbrot.cpp color.cpp)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -fno-rtti -fno-exceptions")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(APPLE)
    target_link_libraries(mandelbrot PRIVATE PNG::PNG)
else()
    target_link_libraries(mandelbrot PRIVATE OpenMP::OpenMP_CXX PNG::PNG)
endif()

