# files
set(DIR_LIB_SRCS 
	pngw_Epart.c
	pngw_LEpart.c
	pngw_Lpart.c
	pngw_Upart.c
	pngw_Utils.c
	pngw.c)

# set dy-library
add_library (core MODULE ${DIR_LIB_SRCS})
add_library (proj::core ALIAS core)
target_include_directories(core PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(core basic)

# profect info
project(MAIN_CORE)

# target
# add_executable(MAIN_CORE main.c)
# target_link_libraries(MAIN_CORE basic)