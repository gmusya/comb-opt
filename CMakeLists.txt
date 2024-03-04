cmake_minimum_required(VERSION 3.22)
project(comb-opt)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ortools CONFIG REQUIRED)

add_subdirectory(tsp)
