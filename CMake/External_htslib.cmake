set(name "htslib")
set(url "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2")
set(dl "${CMAKE_CURRENT_BINARY_DIR}/${name}-dl")
set(src "${CMAKE_CURRENT_BINARY_DIR}/${name}-src")
set(build "${CMAKE_CURRENT_BINARY_DIR}/${name}")
set(install "${CMAKE_CURRENT_BINARY_DIR}/${name}_install")

#get_filename_component(ZLib_cf_LIBRARY_DIR ZLib_cf_LIBRARY DIRECTORY)

set(ZLib_cf_install "${CMAKE_CURRENT_BINARY_DIR}/ZLib_cf_install")
set(ZLib_cf_LIBRARY_DIR ${ZLib_cf_install}/lib)
set(Zlib_cf_INCLUDE ${ZLib_cf_install}/include)

set(HTSLIB_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
set(HTSLIB_CPP_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")

if($ENV{CMAKE} STREQUAL "x86_64-apple-darwin18-cmake")
  set(CONFIG_CROSS_COMPILE_APPLE "--host x86_64-apple-darwin18")
endif()

ExternalProject_Add(
  ${name}_project
  URL ${url}
  DOWNLOAD_DIR ${dl}
  SOURCE_DIR ${src}
  INSTALL_DIR ${install}
  BUILD_IN_SOURCE 1
  #PATCH_COMMAND patch -t -N ${src}/configure.ac ${PROJECT_SOURCE_DIR}/CMake/htslib-macosx.patch
  #CONFIGURE_COMMAND cd ${src} && bash -c "CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} AR=${CMAKE_AR} RANLIB=${CMAKE_RANLIB} rm -v ./configure && autoreconf" && bash -c " AR=${CMAKE_AR} RANLIB=${CMAKE_RANLIB} CC=${CMAKE_C_COMPILER} ./configure --host x86_64-apple-darwin18 --prefix=${install} --without-curses --without-libdeflate --disable-bz2 --disable-lzma --enable-libcurl=no"
  CONFIGURE_COMMAND cd ${src} && bash -c "CFLAGS='-fPIC -I${Zlib_cf_INCLUDE} ${HTSLIB_C_FLAGS}' LDFLAGS='-fPIC -L${ZLib_cf_LIBRARY_DIR}' ./configure ${CONFIG_CROSS_COMPILE_APPLE} CC=${CMAKE_C_COMPILER} --prefix=${install} --without-curses --without-libdeflate --disable-bz2 --disable-lzma --enable-libcurl=no"
#  BUILD_COMMAND cd ${src} && sed -i "s#^CFLAGS.*$#CFLAGS = -fPIC -I${Zlib_cf_INCLUDE} ${HTSLIB_C_FLAGS}#" Makefile && sed -i "s#^LDFLAGS.*$#LDFLAGS = -fPIC -L${ZLib_cf_LIBRARY_DIR}#" Makefile && make
)


# Specify include dir
set(${name}_INCLUDE_DIRS)
list(APPEND ${name}_INCLUDE_DIRS "${src}/htslib-1.9")
set(${name}_INCLUDE_DIR "${${name}_INCLUDE_DIRS}")

# library dir
#if(NOT ${BUILD_SHARED_LIBS})
set(${name}_LIBRARY_PATH ${src}/htslib-1.9/libhts.a)

set(${name}_LIBRARY htslib)
add_library(${${name}_LIBRARY} UNKNOWN IMPORTED)
set_property(TARGET ${${name}_LIBRARY} PROPERTY IMPORTED_LOCATION
                ${${name}_LIBRARY_PATH})

# zlib must be downloaded and compiled before
add_dependencies(${name}_project ${ZLib_cf_LIBRARY})
add_dependencies(${${name}_LIBRARY} ${name}_project)

