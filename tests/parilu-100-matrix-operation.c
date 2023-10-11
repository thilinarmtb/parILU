#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "parilu-impl.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  parilu_set_log_level(3);

  const char *filename = "parilu-100-laplacian-4x4.mtx";
  parilu_matrix *matrix = parilu_matrix_from_file(filename, comm);

  parilu_matrix_operator *matrix_operator =
      parilu_matrix_operator_setup(matrix, comm);

  const uint32_t num_rows = parilu_matrix_num_rows(matrix);
  double *x = (double *)calloc(num_rows, sizeof(double));
  for (uint32_t i = 0; i < num_rows; i++)
    x[i] = 1.0;

  double *y = (double *)calloc(num_rows, sizeof(double));
  parilu_matrix_operator_apply(y, matrix_operator, x);

  // Multiplying Laplacian matrix by \underline{1} should give \underline{0}.
  uint32_t err = 0;
  for (uint32_t i = 0; i < num_rows; i++)
    err |= (fabs(y[i] - 0) > 1e-12);

  free(x), free(y);
  parilu_matrix_operator_free(&matrix_operator);
  parilu_matrix_free(&matrix);

  MPI_Finalize();

  return err;
}
