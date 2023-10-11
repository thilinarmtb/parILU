#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "parilu-impl.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  const char *filename = "parilu-100-mat.mtx";
  uint32_t nnz = 0;
  uint64_t *row = NULL, *col = NULL;
  double *val = NULL;
  parilu_read_matrix(&nnz, &row, &col, &val, filename, comm, 0);

  parilu_matrix *matrix = parilu_matrix_setup(nnz, row, col, val, comm);
  free(row), free(col), free(val);

  parilu_matrix_operator *operator= parilu_matrix_operator_setup(matrix, comm);

  parilu_matrix_operator_free(&operator);
  parilu_matrix_free(&matrix);

  MPI_Finalize();

  return 0;
}
