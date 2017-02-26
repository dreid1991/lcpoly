cmake_minimum_required (VERSION 2.6)
project (lcpoly)
#note that Boost and PYTHON are case-sensitive

file (GLOB SOURCES "src/*cpp")
#REMEMBER OPT FLAGS
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpic -Wno-unused-result -g -O3")


add_executable(lcpoly ${SOURCES})
