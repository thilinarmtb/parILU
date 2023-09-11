#include <parilu.h>
#include <stdio.h>

int main(int argc, char **argv) {
  struct parilu_t *parilu = parilu_init(&argc, &argv);

  for (int i = 0; i < argc; i++)
    printf("%d: %s\n", i, argv[i]);

  parilu_finalize(&parilu);

  return 0;
}
