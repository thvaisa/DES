cmake_minimum_required(VERSION 3.2)

project(DES_PROJECT)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g2")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(DES_INSTALL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(DES_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(DES_INSTALL_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

set(DES_HEADERS_DIR         ${PROJECT_SOURCE_DIR}/src/)
set(DESENGINE_HEADERS_DIR   ${PROJECT_SOURCE_DIR}/src/desengine)
set(GEOMETRY_HEADERS_DIR    ${PROJECT_SOURCE_DIR}/src/desengine/geometry)
set(ORBIT_HEADERS_DIR    ${PROJECT_SOURCE_DIR}/src/desengine/orbit)


include_directories(${DES_INSTALL_INCLUDE_DIR})
include_directories(${DES_HEADERS_DIR})
include_directories(${DESENGINE_HEADERS_DIR})
include_directories(${GEOMETRY_HEADERS_DIR})
include_directories(${ORBIT_HEADERS_DIR})

add_subdirectory(src)
add_subdirectory(test)


