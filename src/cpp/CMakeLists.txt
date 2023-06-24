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


add_subdirectory(../../ext/googletest)
add_subdirectory(tests)

file (GLOB Sources "isotherms/*.cpp")
file (GLOB EmpiricalSources "isotherms/empirical/*.cpp")
file (GLOB EmpiricalFilesSources "isotherms/empirical/files/*.cpp")
file (GLOB OptimizationSources "isotherms/optimization_algorithms/*.cpp")
file (GLOB PotentialFilesSources "isotherms/pta/*.cpp")
file (GLOB PotentialSources "isotherms/pta/*.cpp")
file (GLOB VacancySources "isotherms/vsm/*.cpp")
file (GLOB VacancySourcesActivity "isotherms/vsm/activity_coefficients/*.cpp")
file (GLOB EoSSources "isotherms/equations_of_state/*.cpp")

file (GLOB Headers "isotherms/*.h")
file (GLOB EmpiricalHeaders "isotherms/empirical/*.h")
file (GLOB EmpiricalFilesHeaders "isotherms/empirical/files/*.h")
file (GLOB OptimizationHeaders "isotherms/optimization_algorithms/*.h")
file (GLOB PotentialFilesHeaders "isotherms/pta/*.h")
file (GLOB PotentialHeaders "isotherms/pta/*.h")
file (GLOB VacancyHeaders "isotherms/vsm/*.h")
file (GLOB VacancyHeadersActivity "isotherms/vsm/activity_coefficients/*.h")
file (GLOB EoSHeaders "isotherms/equations_of_state/*.h")




add_library(${This}  SHARED 
    ${Sources} ${OptimizationSources} ${EmpiricalSources} ${EmpiricalFilesSources} ${PotentialSources} ${PotentialFilesSources} ${VacancySources} ${VacancySourcesActivity} ${EoSSources}
    ${Headers} ${OptimizationHeaders} ${EmpiricalHeaders} ${EmpiricalFilesHeaders} ${PotentialHeaders} ${PotentialFilesHeaders} ${VacancyHeaders} ${VacancyHeadersActivity} ${EoSHeaders}
     )