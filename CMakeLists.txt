cmake_minimum_required(VERSION 3.22)
project(comb-opt)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ortools CONFIG REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main ortools::ortools)
