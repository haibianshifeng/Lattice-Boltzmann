add_executable(main main.cu)

target_link_libraries(main PUBLIC ${PROJECT_NAME} ${LIBRARY})

set_target_properties(
        main
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

add_executable(test test.cu)

target_link_libraries(test PUBLIC ${PROJECT_NAME} ${LIBRARY})

set_target_properties(
        test
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)