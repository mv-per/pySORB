cmake_minimum_required(VERSION 3.12)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#  Project
set(This pySorb)

project(${This})

include_directories("${CMAKE_SOURCE_DIR}/pysorb/cpp/isotherms")
include_directories("${CMAKE_SOURCE_DIR}/pysorb/cpp/isotherms/optimization_algorithms")
include_directories("${CMAKE_SOURCE_DIR}/pysorb/cpp/isotherms/equations_of_state")

file (GLOB Sources "pysorb/cpp/isotherms/*.cpp")
file (GLOB EmpiricalSources "pysorb/cpp/isotherms/empirical/*.cpp")
file (GLOB OptimizationSources "pysorb/cpp/isotherms/optimization_algorithms/*.cpp")
file (GLOB PotentialSources "pysorb/cpp/isotherms/pta/*.cpp")
file (GLOB VacancySources "pysorb/cpp/isotherms/vsm/*.cpp")
file (GLOB VacancySourcesActivity "pysorb/cpp/isotherms/vsm/activity_coefficients/*.cpp")
file (GLOB EoSSources "pysorb/cpp/isotherms/equations_of_state/*.cpp")

file (GLOB Headers "pysorb/cpp/isotherms/*.cpp")
file (GLOB EmpiricalHeaders "pysorb/cpp/isotherms/empirical/*.cpp")
file (GLOB OptimizationHeaders "pysorb/cpp/isotherms/optimization_algorithms/*.cpp")
file (GLOB PotentialHeaders "pysorb/cpp/isotherms/pta/*.cpp")
file (GLOB VacancyHeaders "pysorb/cpp/isotherms/vsm/*.cpp")
file (GLOB VacancyHeadersActivity "pysorb/cpp/isotherms/vsm/activity_coefficients/*.cpp")
file (GLOB EoSHeaders "pysorb/cpp/isotherms/equations_of_state/*.cpp")

find_package(Python COMPONENTS Interpreter Development)
add_subdirectory(pybind11)

pybind11_add_module(
    ${This} 
    ${Sources} 
    ${OptimizationSources}
     ${EmpiricalSources} 
     ${PotentialSources} 
     ${VacancySources} 
     ${VacancySourcesActivity} 
     ${EoSSources}
    ${Headers} 
    ${OptimizationHeaders}
     ${EmpiricalHeaders} 
     ${PotentialHeaders} 
     ${VacancyHeaders} 
     ${VacancyHeadersActivity} 
     ${EoSHeaders}
    pysorb/module.cpp)

target_link_libraries(${This} PUBLIC)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(${This}
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})