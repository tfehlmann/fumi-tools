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

ExternalProject_Add(
  ${name}_project
  URL ${url}
  DOWNLOAD_DIR ${dl}
  SOURCE_DIR ${src}
  INSTALL_DIR ${install}
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND cd ${src} && bash -c "CFLAGS='-fPIC -I${Zlib_cf_INCLUDE} ${HTSLIB_C_FLAGS}' LDFLAGS='-fPIC -L${ZLib_cf_LIBRARY_DIR}' ./configure --prefix=${install} --without-curses"
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
add_dependencies(${${name}_LIBRARY} ${name}_project)

