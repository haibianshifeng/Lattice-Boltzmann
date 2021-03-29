aux_source_directory(./app APP)

aux_source_directory(./core BOLTZMANN)

aux_source_directory(./utils UTILS)

add_library(${PROJECT_NAME} SHARED ${APP} ${BOLTZMANN} ${UTILS})

set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)