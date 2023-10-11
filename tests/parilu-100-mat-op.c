#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "parilu-impl.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  parilu_options *options = parilu_default_options();

  const char *filename = "parilu-100-mat.mtx";

  uint32_t nnz = 0;
  uint64_t *row = NULL, *col = NULL;
  double *val = NULL;
  parilu_read_matrix(&nnz, &row, &col, &val, filename, MPI_COMM_WORLD, 0);

  parilu_handle *handle =
      parilu_setup(nnz, row, col, val, options, MPI_COMM_WORLD);

  free(row), free(col), free(val);
  parilu_finalize_options(&options);
  parilu_finalize(&handle);

  MPI_Finalize();

  return 0;
}
