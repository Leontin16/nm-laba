cmake_minimum_required(VERSION 3.10)
project(SplineLab)

set(CMAKE_CXX_STANDARD 17)

add_library(algorithm algorithm.cpp)
target_include_directories(algorithm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(algorithm_test test.cpp)
target_link_libraries(algorithm_test PRIVATE algorithm)

add_executable(accuracy_engine_simple accuracy_engine_simple.cpp)
target_link_libraries(accuracy_engine_simple PRIVATE algorithm)

if(MINGW)
    target_link_libraries(accuracy_engine_simple PRIVATE -static-libstdc++ -static-libgcc)
endif()