set(name "cppformat")
set(url "https://github.com/fmtlib/fmt/releases/download/6.0.0/fmt-6.0.0.zip")
set(dl "${CMAKE_CURRENT_BINARY_DIR}/${name}-dl")
set(src "${CMAKE_CURRENT_BINARY_DIR}/${name}-src")
set(build "${CMAKE_CURRENT_BINARY_DIR}/${name}")
set(install "${CMAKE_CURRENT_BINARY_DIR}/${name}_install")

if(NOT MSVC)
  if(PROCESSOR_COUNT)
    set(parallel_build "-j${PROCESSOR_COUNT}")
  endif()
else()
  set(parallel_build "/MP")
endif()

if(${BUILD_SHARED_LIBS})
  set(POSITION_INDEPENDENT_CODE "-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true")
endif()

ExternalProject_Add(
  ${name}_project
  URL ${url}
  DOWNLOAD_DIR ${dl}
  SOURCE_DIR ${src}
  BINARY_DIR ${build}
  INSTALL_DIR ${install}
  BUILD_COMMAND ""
  CMAKE_ARGS
  "-G${CMAKE_GENERATOR}"
  "-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>"
  "-DCMAKE_INSTALL_LIBDIR:PATH=<INSTALL_DIR>/lib64"
  "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
  '-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} ${CXX_FLAGS}"'
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DFMT_DOC=OFF"
  "-DFMT_TEST=OFF"
  "${POSITION_INDEPENDENT_CODE}"
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
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set(${name}_LIBRARY_PATH ${install}/lib64/libfmtd.a)
else()
  set(${name}_LIBRARY_PATH ${install}/lib64/libfmt.a)
endif()

set(${name}_LIBRARY cppformat)
add_library(${${name}_LIBRARY} UNKNOWN IMPORTED)
set_property(TARGET ${${name}_LIBRARY} PROPERTY IMPORTED_LOCATION
                ${${name}_LIBRARY_PATH})
add_dependencies(${${name}_LIBRARY} ${name}_project)
