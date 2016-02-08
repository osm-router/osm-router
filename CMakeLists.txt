cmake_minimum_required(VERSION 2.8.7)
project(router CXX)
include_directories(include)
ADD_DEFINITIONS(
    -std=c++11 
)
set(SOURCES src/Utils.c++ src/Router.c++)

add_executable(router ${SOURCES})
target_link_libraries(router boost_program_options)
