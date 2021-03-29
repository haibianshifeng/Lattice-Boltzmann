cmake_minimum_required(VERSION 3.17)

project(lattice_boltzmann CUDA)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -g -march=native -fopenmp -std=c++11 -lm")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake_modules/")
set(CMAKE_CUDA_STANDARD 11)

find_package(SFML 2.5 REQUIRED audio graphics window system)
set(LIBRARY sfml-audio sfml-graphics sfml-window sfml-system)
add_subdirectory(src)
target_include_directories(${PROJECT_NAME} PUBLIC include)
add_subdirectory(bin)

add_custom_target(copy-runtime-files ALL COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data)
