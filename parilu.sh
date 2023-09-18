#!/bin/bash

function print_help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  --help Print this help message and exit."
  echo "  --cc <compiler> Set the compiler to use for the build."
  echo "  --build-type <Release|Debug> Build type."
  echo "  --build-dir <build directory> Build directory."
  echo "  --install-prefix <install prefix> Install prefix."
  echo "  --enable-docs Enable building documentation."
  echo "  --enable-asan Enable address sanitizer."
  echo "  --install Install the project."
  echo "  --docs Build documentation with Doxygen+Sphinx."
  echo "  --format Format the source code with clang-format."
  echo "  --format-check Check if source formatting is compliant with clang-format."
  echo "  --tidy Run clang-tidy."
}

# Set default values.
: ${PARILU_CC:=cc}
: ${PARILU_BUILD_TYPE:=Release}
: ${PARILU_INSTALL_PREFIX:=`pwd`/install}
: ${PARILU_BUILD_DIR:=`pwd`/build}
: ${PARILU_ENABLE_DOCS:=OFF}
: ${PARILU_ENABLE_ASAN:=OFF}
: ${PARILU_INSTALL:=NO}
: ${PARILU_DOCS:=NO}
: ${PARILU_FORMAT:=NO}
: ${PARILU_FORMAT_CHECK:=NO}
: ${PARILU_TIDY:=NO}

# Handle command line arguments.
while [[ $# -gt 0 ]]; do
  case $1 in
    --help)
      print_help
      exit 0
      ;;
    --cc)
      PARILU_CC="$2"
      shift
      shift
      ;;
    --build-type)
      PARILU_BUILD_TYPE="$2"
      shift
      shift
      ;;
    --build-dir)
      PARILU_BUILD_DIR="$2"
      shift
      shift
      ;;
    --install-prefix)
      PARILU_INSTALL_PREFIX="$2"
      shift
      shift
      ;;
    --enable-docs)
      PARILU_ENABLE_DOCS="ON"
      shift
      ;;
    --enable-asan)
      PARILU_ENABLE_ASAN="ON"
      shift
      ;;
    --install)
      PARILU_INSTALL="YES"
      shift
      ;;
    --docs)
      PARILU_DOCS="YES"
      shift
      ;;
    --format)
      PARILU_FORMAT="YES"
      shift
      ;;
    --format-check)
      PARILU_FORMAT_CHECK="YES"
      shift
      ;;
    --tidy)
      PARILU_TIDY="YES"
      shift
      ;;
    *)
      echo "Unknown option: $1"
      print_help
      exit 1
      ;;
  esac
done
  
mkdir -p ${PARILU_BUILD_DIR} 2> /dev/null

cmake -DCMAKE_C_COMPILER=${PARILU_CC} \
  -DCMAKE_BUILD_TYPE=${PARILU_BUILD_TYPE} \
  -DCMAKE_INSTALL_PREFIX=${PARILU_INSTALL_PREFIX} \
  -B ${PARILU_BUILD_DIR} \
  -DENABLE_DOCS=${PARILU_ENABLE_DOCS} \
  -DENABLE_ASAN=${PARILU_ENABLE_ASAN} \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -S .
  
if [[ "${PARILU_FORMAT}" == "YES" ]]; then
  cmake --build ${PARILU_BUILD_DIR} --target format -j4
fi

if [[ "${PARILU_FORMAT_CHECK}" == "YES" ]]; then
  cmake --build ${PARILU_BUILD_DIR} --target format-check -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: clang-format check failed."
    exit 1
  fi
fi

if [[ "${PARILU_TIDY}" == "YES" ]]; then
  cmake --build ${PARILU_BUILD_DIR} --target tidy -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: clang-tidy failed."
    exit 1
  fi
fi

if [[ ${PARILU_DOCS} == "YES" ]]; then
  cmake --build ${PARILU_BUILD_DIR} --target docs -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: Building docs with Doxygen+Sphinx failed."
    exit 1
  fi
fi

if [[ "${PARILU_INSTALL}" == "YES" ]]; then
  cmake --build ${PARILU_BUILD_DIR} --target install -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: Installing failed."
    exit 1
  fi
fi
