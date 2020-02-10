cmake_minimum_required(VERSION 3.2)
project(des_tests)

add_subdirectory(lib/googletest)

include_directories(${GEOMETRY_HEADERS_DIR})
include_directories(${DESENGINE_HEADERS_DIR})
include_directories(lib/googletest/googletest/include)

set(SOURCE_FILES main.cpp src/geometry_tests.cpp src/desengine_tests.cpp)


add_executable(des_tests ${SOURCE_FILES})
target_link_libraries(des_tests  desengine geometry gtest)
install(TARGETS des_tests DESTINATION bin)

