add_executable(boltzmann main.cu)

target_link_libraries(boltzmann PUBLIC ${PROJECT_NAME} ${LIBRARY})

set_target_properties(
        boltzmann
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)