cmake_minimum_required(VERSION 3.16)
project(rabbit_qc)

set(CMAKE_CXX_STANDARD 14)
#add_definitions(-DVec512)
add_definitions(-DVec256)
#add_definitions(-DUseAda)
#add_definitions(-DUseLong)
#add_definitions(-DTimer)
#add_definitions(-DKK)
#add_definitions(-D_OPENMP)

aux_source_directory(./src/ DIR_HELLO_SRCS)
add_executable(rabbit_qc ${DIR_HELLO_SRCS})