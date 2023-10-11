#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "parilu.h"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 2 && rank == 0) {
    fprintf(stderr, "Usage: %s <matrix file> [<verbose level>]\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  const char *file = argv[1];
  unsigned verbose = 0;
  if (argc > 2)
    verbose = atoi(argv[2]);

  parilu_options *options = parilu_default_options();
  parilu_set_verbose(options, verbose);

  uint32_t nnz;
  uint64_t *row = NULL, *col = NULL;
  double *val = NULL;
  parilu_read_matrix(&nnz, &row, &col, &val, file, MPI_COMM_WORLD, verbose);

  parilu_handle *ilu =
      parilu_setup(nnz, row, col, val, options, MPI_COMM_WORLD);

  free(row), free(col), free(val);

  parilu_finalize(&ilu);
  parilu_finalize_options(&options);

  MPI_Finalize();

  return 0;
}
