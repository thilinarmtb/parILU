# parILU

## Introduction

parILU - Parallel Incomplete LU factorization.

## Build this project

This project uses conda to manage dependencies (CMake, clang-format, clang-tidy
and other dependencies for documentation). Dependencies can be installed by
executing following commands after installing [conda](https://docs.conda.io/en/latest/miniconda.html).
```sh
conda env create -f environment-dev.yml
conda activate mylib-dev
```

You can format the source code with `clang-format` using the option `--format yes`.
```
./parilu.sh --format yes
```

Then simply run `parilu.sh` script to build and install the library.
```sh
./parilu.sh --docs yes --install yes
```

Use `-h` or `--help` to see all the options supported by `parilu.sh` script.
```
./parilu.sh --help
```
