cmake_minimum_required(VERSION 3.2)
PROJECT(orbit C CXX)

set(SOURCE_FILES
    orbit.hpp
    orbit.cpp
)

add_library(orbit SHARED STATIC ${SOURCE_FILES})

install(TARGETS orbit DESTINATION ${DES_INSTALL_LIB_DIR})
install(FILES orbit.hpp DESTINATION ${DES_INSTALL_INCLUDE_DIR})
