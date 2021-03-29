cmake_minimum_required(VERSION 3.16)

project(lattice_boltzmann CUDA)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake_modules/")
set(CMAKE_CUDA_STANDARD 11)

find_package(SFML 2.5 COMPONENTS graphics window system)
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS})
set(LIBRARY sfml-audio sfml-graphics sfml-window sfml-system ${OPENGL_LIBRARIES})

add_subdirectory(src)
target_include_directories(${PROJECT_NAME} PUBLIC include)
add_subdirectory(bin)
