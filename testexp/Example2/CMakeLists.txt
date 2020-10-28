cmake_minimum_required(VERSION 3.3)

project(UFL C)

option(CXX "enable C++ compilation" ON)
if(CXX)
    enable_language(CXX)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


# list source files here
set(sources Example2.cpp
            Randomclass.cpp)

add_executable(UFL ${sources})

include_directories(${CMAKE_PROJECT_NAME} "/Users/yingqiuzhang/Documents/DSPtest/DSP/src")
include_directories(${CMAKE_PROJECT_NAME} "/Users/yingqiuzhang/Documents/DSPtest/DSP/build/include/coin")
include_directories(${CMAKE_PROJECT_NAME} "/Users/yingqiuzhang/Documents/DSPtest/DSP/testexp/Example2")

target_link_libraries(${CMAKE_PROJECT_NAME} "/Users/yingqiuzhang/Documents/DSPtest/DSP/build/lib/libDsp.dylib")