# files
set(DIR_LIB_SRCS 
	basic_Alloc.c
	basic_Debug.c
	basic_Error.c
	basic_String.c
	utils_Factories.c
	utils_OptionParser.c 
	utils_StringConvert.c
	utils_stringVector.c)

# set dy-library
add_library (basic SHARED ${DIR_LIB_SRCS})
add_library (proj::basic ALIAS basic)
target_include_directories(basic PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(basic)

# profect info
project(MAIN_BASIC)

# target
add_executable(MAIN_BASIC main.c)
target_link_libraries(MAIN_BASIC basic)