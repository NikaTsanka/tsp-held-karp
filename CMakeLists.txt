cmake_minimum_required(VERSION 3.6)
project(tsp_project3)

#set (CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

find_package(X11)
link_libraries(${X11_LIBRARIES})

include_directories(${X11_INCLUDE_DIR})

set(SOURCE_FILES tsp_held_karp.cpp)
add_executable(tsp_project3 ${SOURCE_FILES})