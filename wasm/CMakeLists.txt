cmake_minimum_required(VERSION 3.0)
project(tess_stars2px VERSION 0.1
                  DESCRIPTION "tess_stars2px cmake")

set(CMAKE_BUILD_TYPE DEBUG)

# add the source files needed for main executable
add_executable(tess_stars2px
        src/tess_stars2px.c
        src/mat_ra3.c
        src/vec.c)

# find the math library
find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    # link target with math library
    target_link_libraries(tess_stars2px PUBLIC ${MATH_LIBRARY})
endif()

# define header file path
target_include_directories(tess_stars2px PRIVATE include)

# Make sure the tess_stars2px function is exported along with main for javascript
#  visibility NOTE the leading '_'
target_link_options(tess_stars2px PRIVATE -sEXPORTED_FUNCTIONS=['_main','_tess_stars2px'])

# output wasm and js hopefully
set(CMAKE_EXECUTABLE_SUFFIX ".html")

install(FILES ${CMAKE_BINARY_DIR}/tess_stars2px.js ${CMAKE_BINARY_DIR}/tess_stars2px.wasm DESTINATION ${CMAKE_SOURCE_DIR}/public)
