add_library(algorithm algorithm.cpp)
target_include_directories(algorithm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(algorithm_test test.cpp)
target_link_libraries(algorithm_test PRIVATE algorithm)