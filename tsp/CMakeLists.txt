add_library(tsp
    task.cpp
    task.h
)

target_include_directories(tsp PUBLIC ${CMAKE_SOURCE_DIR})

add_executable(tsp-random random.cpp)
target_link_libraries(tsp-random tsp ortools::ortools)

add_executable(tsp-spanning spanning.cpp)
target_link_libraries(tsp-spanning tsp ortools::ortools)

add_executable(tsp-local local_optimizations.cpp)
target_link_libraries(tsp-local tsp ortools::ortools)

add_executable(tsp-iteratively iteratively.cpp)
target_link_libraries(tsp-iteratively tsp ortools::ortools)

add_executable(tsp-milp-dfj milp_dfj.cpp)
target_link_libraries(tsp-milp-dfj tsp ortools::ortools)

add_executable(tsp-milp-mtz milp_mtz.cpp)
target_link_libraries(tsp-milp-mtz tsp ortools::ortools)

add_executable(tsp-lp-dfj lp_dfj.cpp)
target_link_libraries(tsp-lp-dfj tsp ortools::ortools)

add_executable(tsp-lp-mtz lp_mtz.cpp)
target_link_libraries(tsp-lp-mtz tsp ortools::ortools)
