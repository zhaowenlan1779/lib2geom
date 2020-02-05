cmake_minimum_required(VERSION 3.1)
project(test_dep_2geom CXX C)

option(2GEOM_AS_SUBPROJECT "include 2geom as subproject" OFF)

if (2GEOM_AS_SUBPROJECT)
    message("Using 2geom as subdirectory")
    add_subdirectory("../../" 2geom EXCLUDE_FROM_ALL)
else()
    message("Using installed 2geom")
    find_package(2Geom REQUIRED)
endif()

add_executable(main main.cpp)
target_link_libraries(main PRIVATE 2Geom::2geom)