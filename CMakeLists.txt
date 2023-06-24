cmake_minimum_required(VERSION 3.12)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#  Project
set(This pysorb)

project(${This})

include_directories("${CMAKE_SOURCE_DIR}/src/cpp/isotherms")
include_directories("${CMAKE_SOURCE_DIR}/src/cpp/isotherms/optimization_algorithms")
include_directories("${CMAKE_SOURCE_DIR}/src/cpp/isotherms/equations_of_state")

file (GLOB Sources "src/cpp/isotherms/*.cpp")
file (GLOB EmpiricalSources "src/cpp/isotherms/empirical/*.cpp")
file (GLOB EmpiricalFilesSources "src/cpp/isotherms/empirical/files/*.cpp")
file (GLOB OptimizationSources "src/cpp/isotherms/optimization_algorithms/*.cpp")
file (GLOB PotentialSources "src/cpp/isotherms/pta/*.cpp")
file (GLOB PotentialFilesSources "src/cpp/isotherms/pta/files/*.cpp")
file (GLOB VacancySources "src/cpp/isotherms/vsm/*.cpp")
file (GLOB VacancySourcesActivity "src/cpp/isotherms/vsm/activity_coefficients/*.cpp")
file (GLOB EoSSources "src/cpp/isotherms/equations_of_state/*.cpp")

file (GLOB Headers "src/cpp/isotherms/*.cpp")
file (GLOB EmpiricalHeaders "src/cpp/isotherms/empirical/*.cpp")
file (GLOB EmpiricalFilesHeaders "src/cpp/isotherms/empirical/files/*.cpp")
file (GLOB OptimizationHeaders "src/cpp/isotherms/optimization_algorithms/*.cpp")
file (GLOB PotentialHeaders "src/cpp/isotherms/pta/*.cpp")
file (GLOB PotentialFilesHeaders "src/cpp/isotherms/pta/files/*.cpp")
file (GLOB VacancyHeaders "src/cpp/isotherms/vsm/*.cpp")
file (GLOB VacancyHeadersActivity "src/cpp/isotherms/vsm/activity_coefficients/*.cpp")
file (GLOB EoSHeaders "src/cpp/isotherms/equations_of_state/*.cpp")

find_package(Python COMPONENTS Interpreter Development)
add_subdirectory(ext/pybind11)

pybind11_add_module(
    ${This} 
    ${Sources} 
    ${OptimizationSources}
     ${EmpiricalSources} 
     ${EmpiricalFilesSources} 
     ${PotentialFilesSources} 
     ${PotentialSources} 
     ${VacancySources} 
     ${VacancySourcesActivity} 
     ${EoSSources}
    ${Headers} 
    ${OptimizationHeaders}
     ${EmpiricalHeaders} 
     ${EmpiricalFilesHeaders} 
     ${PotentialHeaders} 
     ${PotentialFilesHeaders} 
     ${VacancyHeaders} 
     ${VacancyHeadersActivity} 
     ${EoSHeaders}
    src/pysorb.cpp)

target_link_libraries(${This} PUBLIC)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(${This}
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})