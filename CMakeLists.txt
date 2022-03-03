# HOW TO USE
# pass the needed compiler if necessary with -DCMAKE_C_COMPILER=/opt/gcc-4.9/bin/gcc -DCMAKE_CXX_COMPILER=/opt/gcc-4.9/bin/g++
# The minimum required version of cmake
cmake_minimum_required(VERSION 2.8)

# The name of the project
project(fumi_tools)

# Use Debug as default build
if(NOT DEFINED CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)

#configure static linking
option(BUILD_SHARED_LIBS "Build a shared (ON) or static (OFF) executable" ON)
option(BUILD_GENERIC "Build so that the code runs on other x86-64 CPUs as well" ON)
option(USE_LIBCPP "Use libc++ instead of libstdc++" OFF)
option(USE_CXXABI "Use cxxabi for clang" OFF)
option(USE_JEMALLOC "Use jemalloc for memory allocation" ON)

if(NOT ${BUILD_SHARED_LIBS})
  #disable -rdynamic
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
  # prefer static libraries
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif(WIN32)
endif()

# Install rpath settings
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# NOT NEEDED RIGHT NOW BECAUSE LIBS ARE STATIC
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1" AND ${BUILD_SHARED_LIBS})
   set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
endif()


# add own modules to the search path
# include add header and add sources
include("${PROJECT_SOURCE_DIR}/CMake/CMakeUtils.cmake")
# include compiler specific additional flags (e.g. {GCC,CLANG,MSVC}_CXX_FLAGS_{DEBUG,RELWITHDEBINFO,RELEASE})
include("${PROJECT_SOURCE_DIR}/CMake/GlobalCompilerFlags.cmake")

set(PRE_CONFIGURE_FILE "${PROJECT_SOURCE_DIR}/include/fumi_tools/version.hpp.in")
set(POST_CONFIGURE_FILE "${PROJECT_SOURCE_DIR}/include/fumi_tools/version.hpp")
set(VERSION_TAG_PREFIX_REGEX "")

if (EXISTS "${PROJECT_SOURCE_DIR}/.git")
    include("${PROJECT_SOURCE_DIR}/CMake/git_watcher.cmake")
elseif(NOT EXISTS "${PROJECT_SOURCE_DIR}/include/fumi_tools/version.hpp")
    message(FATAL_ERROR "The source directory is not a git repository and does not contain a version.hpp file! You have probably downloaded the source code without dependencies. Please download the source code bundled with dependencies.")
endif()


#include("${PROJECT_SOURCE_DIR}/CMake/gitversion/cmake.cmake")

include(ExternalProject)
include("${PROJECT_SOURCE_DIR}/CMake/ExternalZLib.cmake")
include("${PROJECT_SOURCE_DIR}/CMake/External_cppformat.cmake")
include("${PROJECT_SOURCE_DIR}/CMake/External_htslib.cmake")
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${USE_JEMALLOC})
    include("${PROJECT_SOURCE_DIR}/CMake/External_jemalloc.cmake")
else()
    set(JEMALLOC_LIBRARIES "")
endif()


# add the include header tree to the search path for include files
include_directories("${PROJECT_SOURCE_DIR}/include")

# extra library include paths
include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/lib/")
include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/lib/zstr-include")
include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/lib/cxxopts-include")
include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/lib/robin-hood-hashing-include")
include_directories(SYSTEM ${cppformat_INCLUDE_DIR})
include_directories(SYSTEM ${htslib_INCLUDE_DIR})
include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/lib/cpg/include")
include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/lib/string-view-lite/include")
include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/lib/optional-lite/include")

# add code subdirectories
add_subdirectory(include)
add_subdirectory(src)

add_subdirectory(lib/filesystem)

# collect the sources
get_property(SOURCE_FILES GLOBAL PROPERTY SRCS_LIST)
get_property(HEADER_FILES GLOBAL PROPERTY HDRS_LIST)
get_property(MAIN_FILE GLOBAL PROPERTY MAIN_FILE)
get_property(TEST_FILES GLOBAL PROPERTY TESTS_LIST)

#find_package(ZLIB REQUIRED)
set(ZLIB_LIBRARIES ${ZLib_cf_LIBRARY})
set(ZLIB_INCLUDE_DIRS ${ZLib_cf_INCLUDE_DIR})
set(ZLIB_FOUND True)

#include("CMake/TBB/TBBGet.cmake")
#include("CMake/TBB/TBBBuild.cmake")
#tbb_get(TBB_ROOT tbb_root SOURCE_CODE)
#tbb_build(TBB_ROOT ${tbb_root} CONFIG_DIR TBB_DIR)
#find_package(TBB REQUIRED tbb)

if(__COMPILER_GNU AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # GCC, MINGW
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${GCC_CXX_FLAGS_DEBUG}")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -pipe -march=native")
        if(${BUILD_GENERIC})
                set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pipe -march=x86-64 -mtune=generic -msse4.2")
        else()
                set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pipe -march=native")
        endif()

        if(NOT ${BUILD_SHARED_LIBS})
          set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++")
        endif()
        # needed for std::to_string
        if(${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD")
            add_definitions(-D_GLIBCXX_USE_C99 -D_GLIBCXX_USE_C99_MATH -D_GLIBCXX_USE_C99_MATH_TR1)
        endif()
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread")
        if(${USE_LIBCPP})
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
        endif()
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CLANG_CXX_FLAGS_DEBUG}")
				if(${BUILD_GENERIC})
					set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pipe -march=x86-64 -mtune=generic -msse4.2")
				else()
					set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pipe -march=native")
				endif()
        if(NOT ${BUILD_SHARED_LIBS} AND NOT ${APPLE})
          set(CMAKE_EXE_LINKER_FLAGS "-static")
        endif()
endif()

if(MSVC)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${MSVC_CXX_DEBUG_FLAGS}")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MSVC_CXX_RELEASE_FLAGS}")
endif()

if(WIN32)
    add_definitions(-D_WIN32_WINNT=0x0601)
endif()


add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})

add_executable(${PROJECT_NAME}-bin ${POST_CONFIGURE_FILE} ${MAIN_FILE})
add_executable(${PROJECT_NAME}-fix-flags-bin ${POST_CONFIGURE_FILE} src/fix_flags.cpp)
add_executable(${PROJECT_NAME}-demultiplex-bin ${POST_CONFIGURE_FILE} src/demultiplex.cpp src/sample_index_map.cpp)


if (EXISTS "${PROJECT_SOURCE_DIR}/.git")
    add_dependencies(${PROJECT_NAME}-bin check_git_repository)
    add_dependencies(${PROJECT_NAME}-fix-flags-bin check_git_repository)
    add_dependencies(${PROJECT_NAME}-demultiplex-bin check_git_repository)
endif()

set_target_properties(${PROJECT_NAME}-bin
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_dedup)

set_target_properties(${PROJECT_NAME}-fix-flags-bin
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_fix_flags)

set_target_properties(${PROJECT_NAME}-demultiplex-bin
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_demultiplex)


add_dependencies(${PROJECT_NAME} ${cppformat_LIBRARY} ${ZLib_cf_LIBRARY} ${htslib_LIBRARY})

set(COMMON_LIBS ${cppformat_LIBRARY} ${htslib_LIBRARY} ${BZIP2_LIBRARIES} ${LIBLZMA_LIBRARIES} ${CURL_LIBRARIES} ${ZLib_cf_LIBRARY})
# link with libraries
if(NOT WIN32)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${USE_CXXABI})
        set(COMMON_LIBS "c++abi ${COMMON_LIBS}")
    endif()
    target_link_libraries(${PROJECT_NAME} ${COMMON_LIBS} pthread)
    target_link_libraries(${PROJECT_NAME}-bin ${PROJECT_NAME})
    target_link_libraries(${PROJECT_NAME}-fix-flags-bin ${PROJECT_NAME})
    target_link_libraries(${PROJECT_NAME}-demultiplex-bin ${PROJECT_NAME} ${JEMALLOC_LIBRARIES} ghc_filesystem)
else()
    target_link_libraries(${PROJECT_NAME} ${COMMON_LIBS})
    target_link_libraries(${PROJECT_NAME}-bin ${PROJECT_NAME} ws2_32)
    target_link_libraries(${PROJECT_NAME}-demultiplex-bin ${PROJECT_NAME} ws2_32 ${JEMALLOC_LIBRARIES} ghc_filesystem)
endif()

install(TARGETS ${PROJECT_NAME}-bin DESTINATION bin)
install(TARGETS ${PROJECT_NAME}-fix-flags-bin DESTINATION bin)
install(TARGETS ${PROJECT_NAME}-demultiplex-bin DESTINATION bin)

install(PROGRAMS ${PROJECT_SOURCE_DIR}/bin/fumi_tools ${PROJECT_SOURCE_DIR}/bin/fumi_tools_copy_umi DESTINATION bin)
