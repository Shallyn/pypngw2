# version requirement
cmake_minimum_required(VERSION 3.22.2)

# set stdc
set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE INTERNAL "" FORCE)

# project
project(pypngw)

# sub directories
add_subdirectory(libbasic)
add_subdirectory(libcore)
add_subdirectory(libprec)

# target
add_executable(pypngw main.c)
# add_library(pypngw SHARED)
target_link_libraries(pypngw basic)
