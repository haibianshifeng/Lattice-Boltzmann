aux_source_directory(./app APP)

aux_source_directory(./core CORE)

aux_source_directory(./utils UTILS)

add_library(${PROJECT_NAME} SHARED ${APP} ${CORE} ${UTILS})

set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)