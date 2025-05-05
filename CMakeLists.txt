cmake_minimum_required(VERSION 3.16)
project(mpc_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

include_directories(
    ${PROJECT_SOURCE_DIR}/qp_solver/include
    ${PROJECT_SOURCE_DIR}/qp_solver/lib/eigen-3.4.0
)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(qp_module
    src/qp_bindings.cpp
    qp_solver/src/QP.cpp
)

set_target_properties(qp_module PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/src"
)

target_include_directories(qp_module PRIVATE
    ${PROJECT_SOURCE_DIR}/qp_solver/include
    ${PROJECT_SOURCE_DIR}/qp_solver/lib/eigen-3.4.0
)