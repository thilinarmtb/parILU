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

You can format the source code with `clang-format` using the option `--format` if you
make any changes to the source code. Also, you can use `clang-tidy` for static analysis
as well.
```
./parilu.sh --format --tidy
```

Then simply pass build options and targets to be built to `parilu.sh` script:
```sh
./parilu.sh --enable-docs --enable-asan --install
```

Use `--help` to see all the options and targets supported by `parilu.sh` script.
```
./parilu.sh --help
```
