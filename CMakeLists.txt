cmake_minimum_required(VERSION 3.15)
project(StringPool)

set(CMAKE_CXX_STANDARD 14)
#aux_source_directory(./src/ DIR_HELLO_SRCS)
#add_executable(rabbit_qc ${DIR_HELLO_SRCS})

add_executable(StringPool ./src/StringPool.cpp ./src/StringPool.h)