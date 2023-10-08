#include <mpi.h>
#include <stdlib.h>

#include <parilu.h>

static void read_system(uint32_t *const nnz, uint64_t **const row,
                        uint64_t **const col, double **const val,
                        const char *const file) {
  // Only rank 0 will read the file.
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  struct parilu_opts_t *opts = parilu_parse_opts(&argc, &argv);

  free(opts);
  MPI_Finalize();

  return 0;
}
