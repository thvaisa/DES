cmake_minimum_required(VERSION 3.2)
PROJECT(geometry C CXX)

set(SOURCE_FILES
    geometry.hpp
    geometry.cpp
)

add_library(geometry SHARED STATIC ${SOURCE_FILES})

install(TARGETS geometry DESTINATION ${DES_INSTALL_LIB_DIR})
install(FILES geometry.hpp DESTINATION ${DES_INSTALL_INCLUDE_DIR})
