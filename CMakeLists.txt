cmake_minimum_required(VERSION 3.16)

project(lattice_boltzmann CUDA)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake_modules/")
set(CMAKE_CUDA_STANDARD 11)

find_package(SFML 2.5 COMPONENTS graphics window system)
find_package(OpenGL)
include_directories(${OPENGL_INCLUDE_DIRS})

add_subdirectory(src)
target_include_directories(${PROJECT_NAME} PUBLIC include)
add_subdirectory(bin)
