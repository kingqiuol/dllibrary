cmake_minimum_required(VERSION 3.6)

set(CMAKE_CXX_COMPILER D:/vim/mingw/bin/gcc.exe)

set(CMAKE_C_STANDARD 11)

project(dllibrary)


add_subdirectory(image)

add_subdirectory(src)

add_subdirectory(examples)

add_subdirectory(utils)