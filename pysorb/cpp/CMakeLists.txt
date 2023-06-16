cmake_minimum_required(VERSION 3.16)

set(This _pySorb)

project(${This} VERSION 0.1.0 )

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(gtest_force_shared_crt ON)

enable_testing()

# Src
AUX_SOURCE_DIRECTORY(isotherms SRC_FILES)


add_subdirectory(googletest)
add_subdirectory(tests)

file (GLOB Sources "isotherms/*.cpp")
file (GLOB EmpiricalSources "isotherms/empirical/*.cpp")
file (GLOB OptimizationSources "isotherms/optimization_algorithms/*.cpp")
file (GLOB PotentialSources "isotherms/pta/*.cpp")
file (GLOB VacancySources "isotherms/vsm/*.cpp")
file (GLOB VacancySourcesActivity "isotherms/vsm/activity_coefficients/*.cpp")
file (GLOB EoSSources "isotherms/equations_of_state/*.cpp")

file (GLOB Headers "isotherms/*.cpp")
file (GLOB EmpiricalHeaders "isotherms/empirical/*.cpp")
file (GLOB OptimizationHeaders "isotherms/optimization_algorithms/*.cpp")
file (GLOB PotentialHeaders "isotherms/pta/*.cpp")
file (GLOB VacancyHeaders "isotherms/vsm/*.cpp")
file (GLOB VacancyHeadersActivity "isotherms/vsm/activity_coefficients/*.cpp")
file (GLOB EoSHeaders "isotherms/equations_of_state/*.cpp")


add_library(${This}  SHARED 
    ${Sources} ${OptimizationSources} ${EmpiricalSources} ${PotentialSources} ${VacancySources} ${VacancySourcesActivity} ${EoSSources}
    ${Headers} ${OptimizationHeaders} ${EmpiricalHeaders} ${PotentialHeaders} ${VacancyHeaders} ${VacancyHeadersActivity} ${EoSHeaders}
     )