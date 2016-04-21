cmake_minimum_required(VERSION 2.8.7)
project(graph CXX)
#set (CURL_STATICLIB true)
set (TARGET_GRAPH "graph")
set (TARGET_ROUTER "router")

include_directories (include)
ADD_DEFINITIONS(
    -std=c++11 
    #-DCURL_STATICLIB
)
set(SOURCES_GRAPH src/Graph.c++)
set(SOURCES_ROUTER src/Graph.c++)

#add_executable(${TARGET_GRAPH} ${SOURCES_GRAPH})
add_executable(${TARGET_ROUTER} ${SOURCES_ROUTER})

find_path (LIBCURL_INCLUDE_DIR curl/curl.h HINTS ${LIBCURL_DIR}/include)
include_directories(${LIBCURL_INCLUDE_DIR})
#add_dependencies(${TARGET_GRAPH} libcurl)

set (LIBS ${LIBS} boost_program_options boost_filesystem boost_system)
#set (LIBS ${LIBS} curl)
target_link_libraries(${TARGET_ROUTER} ${LIBS})
#target_link_libraries(${TARGET_GRAPH} ${LIBS})
