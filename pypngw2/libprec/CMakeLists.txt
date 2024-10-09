# files
set(DIR_LIB_SRCS 
	pngwprec_Utils.c)

# set dy-library
add_library (prec SHARED ${DIR_LIB_SRCS})
add_library (proj::prec ALIAS prec)
target_include_directories(prec PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(prec basic)

# profect info
project(MAIN_PREC)

# target
# add_executable(MAIN_CORE main.c)
# target_link_libraries(MAIN_CORE basic)