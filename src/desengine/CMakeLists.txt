cmake_minimum_required(VERSION 3.2)
PROJECT(desengine C CXX)

set(SOURCE_FILES
    desengine.hpp
    desengine.cpp
)

add_subdirectory(geometry)
add_subdirectory(orbit)
add_library(desengine SHARED STATIC ${SOURCE_FILES})
target_link_libraries(desengine geometry orbit)


install(TARGETS desengine DESTINATION ${DES_INSTALL_LIB_DIR})
install(FILES desengine.hpp DESTINATION ${DES_INSTALL_INCLUDE_DIR})
