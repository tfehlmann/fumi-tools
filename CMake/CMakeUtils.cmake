# cmake utils
# function to collect all the sources from sub-directories into a single list
function(add_sources)
  get_property(is_defined GLOBAL PROPERTY SRCS_LIST DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY SRCS_LIST
      BRIEF_DOCS "List of source files"
      FULL_DOCS "List of source files to be compiled")
  endif()
  # make absolute paths
  set(SRCS)
  foreach(s IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${s}")
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY SRCS_LIST "${SRCS}")
endfunction(add_sources)

# function to collect all the headers from sub-directories into a single list
function(add_headers)
  get_property(is_defined GLOBAL PROPERTY HDRS_LIST DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY HDRS_LIST
      BRIEF_DOCS "List of header files"
      FULL_DOCS "List of header files")
  endif()
  # make absolute paths
  set(HDRS)
  foreach(s IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND HDRS "${s}")
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY HDRS_LIST "${HDRS}")
endfunction(add_headers)

# function to collect all the sources from sub-directories into a single list
function(add_tests)
  get_property(is_defined GLOBAL PROPERTY TESTS_LIST DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY TESTS_LIST
      BRIEF_DOCS "List of source files"
      FULL_DOCS "List of source files to be compiled")
  endif()
  # make absolute paths
  set(TESTS)
  foreach(s IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND TESTS "${s}")
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY TESTS_LIST "${TESTS}")
endfunction(add_tests)

function(add_main file)
  get_property(is_defined GLOBAL PROPERTY MAIN_FILE DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY MAIN_FILE
      BRIEF_DOCS "File containing the main function."
      FULL_DOCS "File containing the main function")
  endif()
  # make absolute paths
  set(MAIN_FILE)
  if(NOT IS_ABSOLUTE "${file}")
    get_filename_component(file "${file}" ABSOLUTE)
  endif()
  set_property(GLOBAL PROPERTY MAIN_FILE "${file}")
endfunction(add_main)

if(NOT DEFINED PROCESSOR_COUNT)
  # Unknown:
  set(PROCESSOR_COUNT 0)

  # Linux:
  set(cpuinfo_file "/proc/cpuinfo")
  if(EXISTS "${cpuinfo_file}")
    file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
    list(LENGTH procs PROCESSOR_COUNT)
  endif()

  # Mac:
  if(APPLE)
    find_program(cmd_sys_pro "system_profiler")
    if(cmd_sys_pro)
      execute_process(COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info)
      string(REGEX REPLACE "^.*Total Number Of Cores: ([0-9]+).*$" "\\1"
        PROCESSOR_COUNT "${info}")
    endif()
  endif()

  # Windows:
  if(WIN32)
    set(PROCESSOR_COUNT "$ENV{NUMBER_OF_PROCESSORS}")
  endif()
endif()