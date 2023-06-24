cmake_minimum_required(VERSION 3.16)

set(This _pySorbTests)

set(Sources
test_classic_isotherms.cpp
test_vsm.cpp
test_pta.cpp
)

add_executable(${This} ${Sources})
target_link_libraries(${This} PUBLIC
    gtest_main
    _pySorb
)

add_test(
    NAME ${This}
    COMMAND ${This}
)