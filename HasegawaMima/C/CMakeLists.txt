cmake_minimum_required(VERSION 3.12.3)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR})
set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)

project(HasegawaMima)
message(STATUS "Configuring ${PROJECT_NAME}")

add_subdirectory(HM_C)