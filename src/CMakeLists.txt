add_subdirectory(2geom)
add_subdirectory(performance-tests)
if(GTK3_FOUND)
  add_subdirectory(toys)
else()
  message("Not building toys as they require GTK3.")
endif()

add_subdirectory(cython)
add_subdirectory(py2geom)

if (2GEOM_TESTING)
    add_subdirectory(tests)
else()
    message("tests are disabled")
endif()

