cmake_minimum_required(VERSION 2.8.7)
project(graph CXX)
#set (CURL_STATICLIB true)
set (TARGET_GRAPH "graph")
set (TARGET_XML "xml")
set (TARGET_ROUTER "router")
set (TARGET_SP "sp")
set (TARGET_TEST "test")

include_directories (include)
ADD_DEFINITIONS(
    -std=c++11 
    #-DCURL_STATICLIB
)

set(SOURCES_GRAPH src/Graph.cpp)
set(SOURCES_XML src/xml-parser.cpp)
set(SOURCES_ROUTER src/Router.cpp)
set(SOURCES_SP src/stochastic-sp.cpp)
set(SOURCES_TEST src/stochastic-test.cpp)

add_executable(${TARGET_GRAPH} ${SOURCES_GRAPH})
add_executable(${TARGET_XML} ${SOURCES_XML})
add_executable(${TARGET_ROUTER} ${SOURCES_ROUTER})
add_executable(${TARGET_SP} ${SOURCES_SP})
add_executable(${TARGET_TEST} ${SOURCES_TEST})

find_path (LIBCURL_INCLUDE_DIR curl/curl.h HINTS ${LIBCURL_DIR}/include)
include_directories(${LIBCURL_INCLUDE_DIR})
add_dependencies(${TARGET_XML} libcurl)

set (LIBS ${LIBS} boost_program_options boost_filesystem boost_system)
set (LIBS ${LIBS} curl)
target_link_libraries(${TARGET_ROUTER} ${LIBS})
target_link_libraries(${TARGET_GRAPH} ${LIBS})
target_link_libraries(${TARGET_XML} ${LIBS})
target_link_libraries(${TARGET_SP} ${LIBS})
target_link_libraries(${TARGET_TEST} ${LIBS})
