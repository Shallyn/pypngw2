# files
set(DIR_LIB_SRCS 
	interface.c
	pngw3_Utils.c
	pngw3_evolve.c
	pngw3.c
	pngw3_ehfunc.c
	src_hlm/pngw3_l2.c
	src_hlm/pngw3_l3.c
	src_hlm/pngw3_l4.c
	src_hlm/pngw3_l5.c)

# set dy-library
add_library (pngw3 SHARED ${DIR_LIB_SRCS})
add_library (proj::pngw3 ALIAS pngw3)
target_include_directories(pngw3 PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(pngw3 basic pnFourier ${PYPNGW2_MATH_LIB})

# profect info
project(MAIN_PNGW3)

# target
add_executable(MAIN_PNGW3 main.c)
target_link_libraries(MAIN_PNGW3 proj::pngw3 ${PYPNGW2_MATH_LIB})
