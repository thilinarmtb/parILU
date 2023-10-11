cmake_minimum_required(VERSION 3.24)
project(parILU VERSION 0.0.1 DESCRIPTION "" LANGUAGES C)

option(ENABLE_DOCS "Enable Documentation" OFF)
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)

set(PARILU_VERBOSE 0)
set(PARILU_TYPE 1)
set(PARILU_PIVOT 0)
set(PARILU_TOL 1e-6)
set(PARILU_NNZ_PER_ROW 10)

# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
# Use, i.e. don't skip the full RPATH for the build tree.
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# When building, don't use the install RPATH.
set(CMAKE_SKIP_BUILD_RPATH FALSE)
# Set the library location.
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# Add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_COMPILE_WARNING_AS_ERROR ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

# Enable docs if ENABLE_DOCS is set.
if (ENABLE_DOCS)
  add_subdirectory("docs")
endif()

# Enable AddressSanitizer if ENABLE_ASAN is set.
set(PARILU_ADDRESS_SANITIZER FALSE)
if (ENABLE_ASAN)
  include(CheckLinkerFlag)
  check_linker_flag(C -fsanitize=address PARILU_ADDRESS_SANITIZER)
  if (NOT PARILU_ADDRESS_SANITIZER)
    message(FATAL_ERROR "AddressSanitizer is enabled but not supported.")
  endif()
endif()

set(PARILU_SOURCES src/parilu-options.c src/parilu-impl.c src/parilu-setup.c
  src/parilu-mat.c src/parilu-partition.c src/parilu-io.c)
configure_file(${CMAKE_SOURCE_DIR}/include/parilu-defs.h.in
  ${CMAKE_SOURCE_DIR}/include/parilu-defs.h @ONLY)

add_library(parilu SHARED ${PARILU_SOURCES})
set_target_properties(parilu PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1
  PUBLIC_HEADER include/parilu.h)
target_include_directories(parilu PRIVATE include)
target_compile_options(parilu PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>)
target_compile_options(parilu PRIVATE
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-unused-function -Wno-unused-parameter>)
if (PARILU_ADDRESS_SANITIZER)
  target_compile_options(parilu PRIVATE -fsanitize=address)
  target_link_options(parilu PRIVATE -fsanitize=address)
endif()

include(3rd_party/gslib.cmake)

install(TARGETS parilu EXPORT parilu-targets LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
)

macro(parilu_generate_bin_targets target prefix install_dir include_dirs)
  string(REPLACE ${prefix} "" target_exe ${target})
  string(REPLACE ".c" "" target_exe ${target_exe})
  set(target_exe ${target_exe})

  add_executable(${target_exe} ${target})
  target_link_libraries(${target_exe} PRIVATE parilu)
  target_include_directories(${target_exe} PRIVATE include ${include_dirs})
  target_compile_options(${target_exe} PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>)
  target_compile_options(${target_exe} PRIVATE
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-unused-function -Wno-unused-parameter>)
  if (PARILU_ADDRESS_SANITIZER)
    target_compile_options(${target_exe} PRIVATE -fsanitize=address)
    target_link_options(${target_exe} PRIVATE -fsanitize=address)
  endif()

  install(TARGETS ${target_exe} RUNTIME DESTINATION ${install_dir})
endmacro(parilu_generate_bin_targets)

# Add examples to the build tree.
file(GLOB PARILU_EXAMPLES ${CMAKE_SOURCE_DIR}/examples/*.c)
foreach(example ${PARILU_EXAMPLES})
  parilu_generate_bin_targets(${example} "${CMAKE_SOURCE_DIR}/examples/"
    "examples" "")
endforeach()

# Add tests to the build tree.
file(GLOB PARILU_TESTS ${CMAKE_SOURCE_DIR}/tests/*.c)
foreach(test ${PARILU_TESTS})
  parilu_generate_bin_targets(${test} "${CMAKE_SOURCE_DIR}/tests/"
    "tests" ${GSLIB_INSTALL_DIR}/include)
endforeach()

# Add clang-format as a custom target if available.
find_program(CLANG_FORMAT NAMES clang-format)
if (CLANG_FORMAT)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT} -i ${CMAKE_CURRENT_SOURCE_DIR}/*/*.[ch]
    COMMENT "Running clang-format ...")
  add_custom_target(format-check
    COMMAND ${CLANG_FORMAT} --dry-run -Werror -i ${CMAKE_CURRENT_SOURCE_DIR}/*/*.[ch]
    COMMENT "Running clang-format check ...")
endif()

# Configure .clang-tidy.in to generate .clang-tidy.
configure_file(${CMAKE_SOURCE_DIR}/.clang-tidy.in ${CMAKE_SOURCE_DIR}/.clang-tidy @ONLY)
# Add clang-tidy as a custom target if available.
find_program(CLANG_TIDY NAMES clang-tidy)
if (CLANG_TIDY)
  add_custom_target(tidy
    COMMAND ${CLANG_TIDY} -p ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*/*.[ch]
    COMMENT "Running clang-tidy ...")
endif()
