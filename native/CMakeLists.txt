cmake_minimum_required(VERSION 3.20)
project(molir_native LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Placeholder for native chemical graph isomorphism matcher (VF3 / RDKit C++)
add_library(molir_graph_matcher STATIC
    graph_matcher/matcher.cpp
)

target_include_directories(molir_graph_matcher PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/graph_matcher
)
