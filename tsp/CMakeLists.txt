add_library(tsp
    task.cpp
    task.h
)

target_include_directories(tsp PUBLIC ${CMAKE_SOURCE_DIR})

add_executable(tsp-random random.cpp)
target_link_libraries(tsp-random tsp ortools::ortools)
