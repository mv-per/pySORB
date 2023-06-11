cmake_minimum_required(VERSION 3.16)

set(This cClassicIsotherms)

project(${This} VERSION 0.1.0 )

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(gtest_force_shared_crt ON)

enable_testing()

# Src
AUX_SOURCE_DIRECTORY(src SRC_FILES)


file (GLOB Sources "src/*.cpp")
file (GLOB Headers "src/*.h")


add_library(${This} SHARED ${Sources} ${Headers})