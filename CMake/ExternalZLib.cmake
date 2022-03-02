set(name "ZLib_cf")
set(url "https://github.com/cloudflare/zlib/archive/ab50ae879a9e9b2d5874b7da17ec6e8b3eba2b85.zip")
set(fname "ZLib_cf.zip")
set(dl "${CMAKE_CURRENT_BINARY_DIR}/${name}-dl")
set(src "${CMAKE_CURRENT_BINARY_DIR}/${name}-src")
set(build "${CMAKE_CURRENT_BINARY_DIR}/${name}")
set(install "${CMAKE_CURRENT_BINARY_DIR}/${name}_install")

if(__COMPILER_GNU) # GCC, MINGW, CLANG
  if(${BUILD_GENERIC})
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -pipe -march=x86-64 -mtune=generic -msse4.2")
  else()
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -pipe -march=native")
  endif()
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -pipe -march=native")
endif()


if(NOT MSVC)
  if(PROCESSOR_COUNT)
    set(parallel_build "-j${PROCESSOR_COUNT}")
  endif()
else()
  set(parallel_build "/MP")
endif()

set(ZLIB_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")

ExternalProject_Add(
  ${name}_project
  URL ${url}
  DOWNLOAD_DIR ${dl}
  SOURCE_DIR ${src}
  BINARY_DIR ${build}
  INSTALL_DIR ${install}
  BUILD_COMMAND ""
  PATCH_COMMAND ""
  CMAKE_ARGS
  "-G${CMAKE_GENERATOR}"
  "-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>"
  "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
  "-DCMAKE_C_FLAGS=${ZLIB_C_FLAGS}"
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DBUILD_SHARED_LIBS=OFF"
  "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}"
  INSTALL_COMMAND
    "${CMAKE_COMMAND}"
    --build .
    --target install
    --config ${CMAKE_BUILD_TYPE}
    -- "${parallel_build}"
)


# Specify include dir
set(${name}_INCLUDE_DIRS)
list(APPEND ${name}_INCLUDE_DIRS "${install}/include")
set(${name}_INCLUDE_DIR "${${name}_INCLUDE_DIRS}")

# library dir
set(${name}_LIBRARY_PATH ${install}/lib/libz.a)

set(${name}_LIBRARY ZLib_cf)
add_library(${${name}_LIBRARY} UNKNOWN IMPORTED)
set_property(TARGET ${${name}_LIBRARY} PROPERTY IMPORTED_LOCATION
                ${${name}_LIBRARY_PATH})
add_dependencies(${${name}_LIBRARY} ${name}_project)

