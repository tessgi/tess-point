cmake_minimum_required(VERSION 3.0)
project(starspx6 VERSION 0.1
                  DESCRIPTION "tess_stars2px cmake"
                  LANGUAGES C)

# Set c Compiler
set(CMAKE_C_COMPILER clang)

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
