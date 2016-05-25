cmake_minimum_required(VERSION 3.3)
project(SPHSm6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp Particile.cpp Particile.h Flow.cpp Flow.h Help.cpp Help.h)
add_executable(SPHSm6 ${SOURCE_FILES})