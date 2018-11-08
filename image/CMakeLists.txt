cmake_minimum_required(VERSION 2.8)
set(CMAKE_C_STANDARD 11)

project(image)

include_directories("..")


set(HEARS
        image.h
        matrix.h
        spatial_filter.h
        panorama_image.h
        stb_image.h
        stb_image_write.h)

set(SOURCE_DIR
        image.c
        matrix.c
        panorama_image.c
        spatial_filter.c)

#include_directories(${PROJECT_SOURCE_DIR}/)
#add_executable(image  ${SOURCE_DIR})

add_library(image ${HEARS} ${SOURCE_DIR})