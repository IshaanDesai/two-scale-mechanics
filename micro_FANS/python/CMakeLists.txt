cmake_minimum_required(VERSION 3.5)

# ##############################################################################
# GENERAL SETTINGS
# ##############################################################################

project(MicroFANS
        VERSION 0.0.1
        LANGUAGES C CXX
)

include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-fmax-errors=1")

FetchContent_Declare(
  pybind11 GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.12.0
)
FetchContent_MakeAvailable(pybind11)

find_package(FANS REQUIRED)

pybind11_add_module(MicroFANS micro.hpp micro.cpp)
target_link_libraries(MicroFANS PRIVATE FANS::FANS)

pybind11_add_module(MicroFANSTHERMAL micro_thermal.hpp micro_thermal.cpp)
target_link_libraries(MicroFANSTHERMAL PRIVATE FANS::FANS)


add_custom_command(
        TARGET MicroFANS
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
        $<TARGET_FILE:MicroFANS>
        ${CMAKE_CURRENT_SOURCE_DIR}/../test/$<TARGET_FILE_NAME:MicroFANS>
        COMMENT "Create a symlink for FANS python bindings to ${CMAKE_CURRENT_SOURCE_DIR}/../test/"
)

add_custom_command(
        TARGET MicroFANSTHERMAL
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
        $<TARGET_FILE:MicroFANSTHERMAL>
        ${CMAKE_CURRENT_SOURCE_DIR}/../test/$<TARGET_FILE_NAME:MicroFANSTHERMAL>
        COMMENT "Create a symlink for FANS python bindings to ${CMAKE_CURRENT_SOURCE_DIR}/../test/"
)

option(USE_MPI "Use MPI" OFF)

if (USE_MPI)
        add_definitions(-DUSE_MPI)
endif()