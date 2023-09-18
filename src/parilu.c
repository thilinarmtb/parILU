#include "parilu-impl.h"
#include <getopt.h>

enum parilu_input_t {
  PARILU_INPUT_VERBOSE = 0,
  PARILU_INPUT_TYPE = 10,
  PARILU_INPUT_TOL = 20,
  PARILU_INPUT_NNZ_PER_ROW = 30,
  PARILU_INPUT_HELP = 99
};

#define PARILU_VERBOSE 0
#define PARILU_TYPE 1
#define PARILU_TOL 1e-6
#define PARILU_NNZ_PER_ROW 10

static void parilu_print_help(const char *name) {
  printf("Usage: %s [OPTIONS]\n", name);
  printf("Options:\n");
  printf("  --parilu-verbose=<verbose level>, Verbose level: 0, 1, "
         "2, ... (%d, Default)\n",
         PARILU_VERBOSE);
  printf("  --parilu-type=<type>, ILU type: 0 (ILU0), 1 (ILUt) (%d, Default)\n",
         PARILU_TYPE);
  printf("  --parilu-tol=<tolerance>, Tolerance (%f, Default)\n", PARILU_TOL);
  printf("  --nnz-per-row=<nnz>, Number of non-zeros per row: 1, 2, ... (%d, "
         "Default)\n",
         PARILU_NNZ_PER_ROW);
  printf("  --parilu-help, Prints this help message and exit.\n");
}

static inline int parilu_parse_number(void **out, const char *const arg,
                                      const uint type) {
  union {
    uint u;
    double d;
  } u;

  char *end = NULL;
  switch (type) {
  case 0:
    u.u = (uint)strtol(arg, &end, 10);
    memcpy(*out, &u.u, sizeof(uint));
    break;
  case 1:
    u.d = strtod(arg, &end);
    memcpy(*out, &u.d, sizeof(double));
    break;
  default:
    break;
  }

  if (!end || *end != '\0' || arg == end)
    return 1;
  return 0;
}

static void parilu_parse_opts(struct parilu_opts_t *parilu, int *argc,
                              char ***argv_) {
  parilu->verbose = PARILU_VERBOSE;
  parilu->type = PARILU_TYPE;
  parilu->tol = PARILU_TOL;
  parilu->nnz_per_row = PARILU_NNZ_PER_ROW;

  static struct option long_options[] = {
      {"parilu-verbose", optional_argument, 0, PARILU_INPUT_VERBOSE},
      {"parilu-type", optional_argument, 0, PARILU_INPUT_TYPE},
      {"parilu-help", no_argument, 0, PARILU_INPUT_HELP},
      {0, 0, 0, 0}};

  char **argv = *argv_;
  void *out[] = {(void *)&parilu->verbose, (void *)&parilu->type,
                 (void *)&parilu->tol, (void *)&parilu->nnz_per_row};
  for (;;) {
    int opt = getopt_long(*argc, argv, "", long_options, NULL);
    if (opt == -1)
      break;

    switch (opt) {
    case PARILU_INPUT_VERBOSE:
      if (parilu_parse_number(&out[0], optarg, 0))
        parilu->verbose = PARILU_VERBOSE;
      break;
    case PARILU_INPUT_TYPE:
      if (parilu_parse_number(&out[1], optarg, 0))
        parilu->type = PARILU_TYPE;
      break;
    case PARILU_INPUT_TOL:
      if (parilu_parse_number(&out[2], optarg, 1))
        parilu->tol = PARILU_TOL;
      break;
    case PARILU_INPUT_NNZ_PER_ROW:
      if (parilu_parse_number(&out[3], optarg, 0))
        parilu->nnz_per_row = PARILU_NNZ_PER_ROW;
      break;
    case PARILU_INPUT_HELP:
      parilu_print_help(argv[0]);
      exit(EXIT_SUCCESS);
      break;
    default:
      parilu_print_help(argv[0]);
      exit(EXIT_FAILURE);
      break;
    }
  }

  // Remove parsed arguments from argv. We just need to update the pointers
  // since command line arguments are not transient and available until the
  // end of the program.
  for (int i = optind; i < *argc; i++)
    argv[i - optind] = argv[i];
  *argc -= optind;
}

/**
 * @ingroup parilu_user_api_functions
 *
 * @brief Initialize parilu_opts_t object from command line input. Returns a
 * pointer to a newly allocated struct parilu_opts_t.
 *
 * @param argc Pointer to the number of commnad line arguments.
 * @param argv Pointer to the array of command line arguments.
 *
 * @return struct parilu_opts_t*
 */
struct parilu_opts_t *parilu_init(int *argc, char **argv[]) {
  struct parilu_opts_t *opts = parilu_calloc(struct parilu_opts_t, 1);

  parilu_parse_opts(opts, argc, argv);

  return opts;
}

/**
 * @ingroup parilu_user_api_functions
 *
 * @brief Finalize parilu. Frees the memory allocated for the struct parilu_t
 * returned by parilu_init() and sets the pointer to NULL.
 *
 * @param parilu Pointer to the struct parilu_t* to be freed.
 */
void parilu_finalize(struct parilu_t **parilu) { parilu_free(parilu); }

#undef PARILU_VERBOSE
#undef PARILU_TYPE
#undef PARILU_TOL
#undef PARILU_NNZ_PER_ROW
