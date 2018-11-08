project(utils)

include_directories("..")

set(HEADER
        utils.h
        blas.h
        im2col.h
        gemm.h
        list.h)

set(SOURCE_DIR
        utils.c
        blas.c
        im2col.c
        gemm.c
        list.c)

add_library(utils ${HEADER} ${SOURCE_DIR})