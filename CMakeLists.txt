cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(CMakeSuperLUAll VERSION 5.2.2 LANGUAGES C)

include(${PROJECT_SOURCE_DIR}/cmake/setup_superlu.cmake)

add_subdirectory(superlu)