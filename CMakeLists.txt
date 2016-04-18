cmake_minimum_required(VERSION 2.8.7)
project(osm-router CXX)
include_directories(include)
ADD_DEFINITIONS(
    -std=c++11 
)
set(SOURCES-ROUTER src/Router.c++)
set(SOURCES-GRAPH src/Graph.c++)

add_executable(router ${SOURCES-ROUTER})
add_executable(graph ${SOURCES-GRAPH})
target_link_libraries(router boost_program_options boost_filesystem boost_system)
target_link_libraries(graph boost_program_options boost_filesystem boost_system)
