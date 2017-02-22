cmake_minimum_required (VERSION 2.6)
project (lcpoly)
#note that Boost and PYTHON are case-sensitive

file (GLOB SOURCES "src/*cpp")
#REMEMBER OPT FLAGS
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpic ")

set (NETWORKS_LIB_NAME Sim)

#add_executable(main_exec main.cpp)
#target_link_libraries (main_exec ${NETWORKS_LIB_NAME})

#build python shared library

add_executable(lcpoly ${SOURCES})
