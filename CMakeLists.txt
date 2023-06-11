cmake_minimum_required(VERSION 3.12)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


#  Project
set(This ClassicIsotherms)

project(${This})

include_directories("${CMAKE_SOURCE_DIR}/isotherms/cpp/src")


file (GLOB Sources "isotherms/cpp/src/*.cpp")
file (GLOB Headers "isotherms/cpp/src/*.h")

find_package(Python COMPONENTS Interpreter Development)
add_subdirectory(pybind11)
pybind11_add_module(
    ${This} 
    ${Sources} 
    ${Headers} 
    isotherms/module.cpp)

target_link_libraries(${This} PUBLIC)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(${This}
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})