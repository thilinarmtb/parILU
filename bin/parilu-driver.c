#include <parilu.h>
#include <stdio.h>

int main(int argc, char **argv) {
  struct parilu_opts_t *opts = parilu_init(&argc, &argv);
  if (opts)
    free(opts);
  return 0;
}
