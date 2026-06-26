# version requirement
cmake_minimum_required(VERSION 3.22.2)

# project
project(pypngw)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)

if(APPLE)
    set(PYPNGW2_NATIVE_RPATH
        "@loader_path"
        "@loader_path/../libbasic"
        "@loader_path/../libcore"
        "@loader_path/../libpnfourier"
        "@loader_path/../libpngw3")
else()
    set(PYPNGW2_NATIVE_RPATH
        "$ORIGIN"
        "$ORIGIN/../libbasic"
        "$ORIGIN/../libcore"
        "$ORIGIN/../libpnfourier"
        "$ORIGIN/../libpngw3")
endif()

set(CMAKE_BUILD_RPATH "${PYPNGW2_NATIVE_RPATH}")
set(CMAKE_INSTALL_RPATH "${PYPNGW2_NATIVE_RPATH}")

# sub directories
add_subdirectory(libbasic)
add_subdirectory(libcore)
add_subdirectory(libprec)
add_subdirectory(libpnfourier)
add_subdirectory(libpngw3)

# target
add_executable(pypngw main.c)
# add_library(pypngw SHARED)
target_link_libraries(pypngw basic)
