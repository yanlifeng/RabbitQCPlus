cmake_minimum_required(VERSION 3.15)
project(StringPool)

set(CMAKE_CXX_STANDARD 14)

add_executable(StringPool
        sequence.cpp sequence.h StringPool.cpp StringPool.h)
