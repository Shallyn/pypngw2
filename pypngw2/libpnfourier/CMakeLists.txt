# files
set(DIR_LIB_SRCS 
	pnFourier_utils.c
	pnFourier_Jpqa0.c
	pnFourier_Jpqab.c
	pnFourier_Kpqab.c
	pnFourier.c)

# set dy-library
add_library (pnFourier SHARED ${DIR_LIB_SRCS})
add_library (proj::pnFourier ALIAS pnFourier)
target_include_directories(pnFourier PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(pnFourier basic gsl fftw3 ${PYPNGW2_MATH_LIB})

# profect info
# project(MAIN_BASIC)

# target
# add_executable(MAIN_BASIC main.c)
# target_link_libraries(MAIN_BASIC basic)
