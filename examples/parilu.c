#include <parilu.h>
#include <stdio.h>

static void read_system(uint32_t *nnz, uint64_t **row, uint64_t **col,
                        double **val) {}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  struct parilu_opts_t *opts = parilu_parse_opts(&argc, &argv);

  free(opts);
  MPI_Finalize();

  return 0;
}
