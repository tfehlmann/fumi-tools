set(name "jemalloc")
set(url "https://github.com/jemalloc/jemalloc/releases/download/5.2.1/jemalloc-5.2.1.tar.bz2")
set(dl "${CMAKE_CURRENT_BINARY_DIR}/${name}-dl")
set(src "${CMAKE_CURRENT_BINARY_DIR}/${name}-src")
set(build "${CMAKE_CURRENT_BINARY_DIR}/${name}")
set(install "${CMAKE_CURRENT_BINARY_DIR}/${name}_install")

ExternalProject_Add(
  ${name}_project
  URL ${url}
  DOWNLOAD_DIR ${dl}
  SOURCE_DIR ${src}
  INSTALL_DIR ${install}
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND cd ${src} && bash -c "./configure --prefix=${install}"
)

# library dir
set(${name}_LIBRARY_PATH ${install}/lib/libjemalloc_pic.a)

set(${name}_LIBRARY jemalloc)
add_library(${${name}_LIBRARY} UNKNOWN IMPORTED)
set_property(TARGET ${${name}_LIBRARY} PROPERTY IMPORTED_LOCATION
                ${${name}_LIBRARY_PATH})

add_dependencies(${${name}_LIBRARY} ${name}_project)

