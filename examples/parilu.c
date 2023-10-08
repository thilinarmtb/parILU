#include <parilu.h>
#include <stdio.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  struct parilu_opts_t *opts = parilu_parse_opts(&argc, &argv);

  if (opts)
    free(opts);

  MPI_Finalize();

  return 0;
}
