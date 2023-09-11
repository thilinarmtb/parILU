#!/bin/bash

function print_help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  -c|--cc <compiler> Set the compiler to use for the build."
  echo "  -t|--type <Release|Debug> Build type."
  echo "  -p|--prefix <install prefix> Install prefix."
  echo "  -b|--build-dir <build directory> Build directory."
  echo "  -d|--docs <yes|no> Enable or disable building documentation."
  echo "  -a|--asan <yes|no> Enable or disable address sanitizer."
  echo "  -h|--help Print this help message and exit."
  echo "  --install <yes|no> Install the project."
  echo "  --format <yes|no> Format the source code with clang-format."
  echo "  --format-check <yes|no> Check if source formatting is compliant with clang-format."
  echo "  --tidy <yes|no> Run clang-tidy."
}

# Set default values.
: ${PARILU_CC:=cc}
: ${PARILU_BUILD_TYPE:=Release}
: ${PARILU_INSTALL_PREFIX:=`pwd`/install}
: ${PARILU_BUILD_DIR:=`pwd`/build}
: ${PARILU_ENABLE_DOCS:=no}
: ${PARILU_ENABLE_ASAN:=no}
: ${PARILU_INSTALL:=no}
: ${PARILU_FORMAT:=no}
: ${PARILU_FORMAT_CHECK:=no}
: ${PARILU_TIDY:=no}

# Handle command line arguments.
while [[ $# -gt 0 ]]; do
  case $1 in
    -c|--cc)
      PARILU_CC="$2"
      shift
      shift
      ;;
    -t|--type)
      PARILU_BUILD_TYPE="$2"
      shift
      shift
      ;;
    -p|--prefix)
      PARILU_INSTALL_PREFIX="$2"
      shift
      shift
      ;;
    -b|--build-dir)
      PARILU_BUILD_DIR="$2"
      shift
      shift
      ;;
    -d|--docs)
      PARILU_ENABLE_DOCS="$2"
      shift
      shift
      ;;
    -a|--asan)
      PARILU_ENABLE_ASAN="$2"
      shift
      shift
      ;;
    -h|--help)
      print_help
      exit 0
      ;;
    --install)
      PARILU_INSTALL="$2"
      shift
      shift
      ;;
    --format)
      PARILU_FORMAT="$2"
      shift
      shift
      ;;
    --format-check)
      PARILU_FORMAT_CHECK="$2"
      shift
      shift
      ;;
    --tidy)
      PARILU_TIDY="$2"
      shift
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
  
if [[ "${PARILU_FORMAT}" == "yes" ]]; then
  cmake --build ${PARILU_BUILD_DIR} --target format -j4
fi

if [[ "${PARILU_FORMAT_CHECK}" == "yes" ]]; then
  cmake --build ${PARILU_BUILD_DIR} --target format-check -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: clang-format check failed."
    exit 1
  fi
fi

if [[ "${PARILU_TIDY}" == "yes" ]]; then
  cmake --build ${PARILU_BUILD_DIR} --target tidy -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: clang-tidy failed."
    exit 1
  fi
fi

if [[ ${PARILU_ENABLE_DOCS} == "yes" ]]; then
  cmake --build ${PARILU_BUILD_DIR} --target Sphinx -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: Building docs with Sphinx failed."
    exit 1
  fi
fi

if [[ "${PARILU_INSTALL}" == "yes" ]]; then
  cmake --build ${PARILU_BUILD_DIR} --target install -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: Installing failed."
    exit 1
  fi
fi
