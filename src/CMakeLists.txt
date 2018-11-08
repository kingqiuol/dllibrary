project(src)

include_directories("..")

set(HEADER_DIR
        convolutional_layer.h
        layer.h
        batchnorm_layer.h
        network.h
        activations.h)

set(SOURCE_DIR
        convolutional_layer.c
        layer.c
        batchnorm_layer.c
        network.c
        activations.c)

add_library(src ${HEADER_DIR} ${SOURCE_DIR})
target_link_libraries(src image utils)
