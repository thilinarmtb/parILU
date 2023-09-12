include(ExternalProject)

set(GSLIB_SOURCE_DIR ${CMAKE_BINARY_DIR}/gslib)
set(GSLIB_INSTALL_DIR ${GSLIB_SOURCE_DIR}/build)
set(GSLIB_CFLAGS -fPIC ${CMAKE_C_FLAGS})

ExternalProject_Add(gslib
  GIT_REPOSITORY    https://github.com/Nek5000/gslib.git
  GIT_TAG           v1.0.8
  SOURCE_DIR        ${GSLIB_SOURCE_DIR}
  INSTALL_DIR       ${GSLIB_INSTALL_DIR}
  BUILD_IN_SOURCE   TRUE
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     CC=${CMAKE_C_COMPILER} LIBNAME=gs CFLAGS=${GSLIB_CFLAGS} make -C <SOURCE_DIR>
)

add_dependencies(parilu gslib)
target_link_libraries(parilu ${GSLIB_INSTALL_DIR}/lib/libgs.a)
target_include_directories(parilu PRIVATE ${GSLIB_INSTALL_DIR}/include)
target_include_directories(parilu-driver PRIVATE ${GSLIB_INSTALL_DIR}/include)
