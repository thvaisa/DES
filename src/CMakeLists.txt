cmake_minimum_required(VERSION 3.2)
PROJECT(des)
enable_language(C CXX)

add_subdirectory(desengine)
set(SOURCE_FILES main.cpp)

add_executable(des ${SOURCE_FILES})

target_link_libraries(des desengine)
install(TARGETS des DESTINATION ${DES_INSTALL_BIN_DIR})
