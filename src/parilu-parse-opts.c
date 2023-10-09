#include <errno.h>
#include <getopt.h>

#include "parilu-impl.h"

enum parilu_input_t {
  PARILU_INPUT_VERBOSE = 0,
  PARILU_INPUT_TYPE = 10,
  PARILU_INPUT_PIVOT = 20,
  PARILU_INPUT_TOL = 40,
  PARILU_INPUT_NNZ_PER_ROW = 50,
  PARILU_INPUT_NULL_SPACE = 60,
  PARILU_INPUT_FILE = 70,
  PARILU_INPUT_HELP = 99
};

static void parilu_print_help(const char *name) {
  printf("Usage: %s [OPTIONS]\n", name);
  printf("Options:\n");
  printf("  --parilu-verbose=<verbose level>, Verbose level: 0, 1, "
         "2, ... (%d, Default).\n",
         PARILU_VERBOSE);
  printf(
      "  --parilu-type=<type>, ILU type: 0 (ILU0), 1 (ILUt) (%d, Default).\n",
      PARILU_TYPE);
  printf("  --parilu-pivot=<pivot>, Pivot: 0 (No pivot), 1 (Partial pivot) "
         "(%d, Default).\n",
         PARILU_PIVOT);
  printf("  --parilu-tol=<tolerance>, Tolerance (%f, Default).\n", PARILU_TOL);
  printf("  --nnz-per-row=<nnz>, Number of non-zeros per row: 1, 2, ... (%d, "
         "Default).\n",
         PARILU_NNZ_PER_ROW);
  printf("  --parilu-file <file>, Input file name with matrix.\n");
  printf("  --parilu-help, Prints this help message and exit.\n");
}

union number {
  int64_t i;
  double d;
};

static void parilu_parse_number(void *out, const char *const arg,
                                const uint type, union number default_value) {
  union number u;
  char *end = NULL;
  size_t unit_size;
  errno = 0;
  switch (type) {
  case 0:
    // NOLINTNEXTLINE
    u.i = strtol(arg, &end, 10);
    unit_size = sizeof(uint);
    memcpy(out, &u.i, unit_size);
    break;
  case 1:
    u.d = strtod(arg, &end);
    unit_size = sizeof(double);
    memcpy(out, &u.d, unit_size);
    break;
  default:
    break;
  }

  if (errno != 0 || *end != '\0' || arg == end)
    memcpy(out, &default_value, unit_size);
}

static void parilu_parse_opts_aux(struct parilu_opts_t *parilu, const int *argc,
                                  char ***argv_) {
  parilu->verbose = PARILU_VERBOSE;
  parilu->type = PARILU_TYPE;
  parilu->pivot = PARILU_PIVOT;
  parilu->tol = PARILU_TOL;
  parilu->nnz_per_row = PARILU_NNZ_PER_ROW;

  static struct option long_options[] = {
      {"parilu-verbose", optional_argument, 0, PARILU_INPUT_VERBOSE},
      {"parilu-type", optional_argument, 0, PARILU_INPUT_TYPE},
      {"parilu-pivot", optional_argument, 0, PARILU_INPUT_PIVOT},
      {"parilu-tol", optional_argument, 0, PARILU_INPUT_TOL},
      {"parilu-nnz-per-row", optional_argument, 0, PARILU_INPUT_NNZ_PER_ROW},
      {"parilu-null-space", required_argument, 0, PARILU_INPUT_NULL_SPACE},
      {"parilu-file", required_argument, 0, PARILU_INPUT_FILE},
      {"parilu-help", no_argument, 0, PARILU_INPUT_HELP},
      {0, 0, 0, 0}};

  union number default_value;
  char **argv = *argv_;
  for (;;) {
    int opt = getopt_long(*argc, argv, "", long_options, NULL);
    if (opt == -1)
      break;

    switch (opt) {
    case PARILU_INPUT_VERBOSE:
      default_value.i = PARILU_VERBOSE;
      parilu_parse_number(&parilu->verbose, optarg, 0, default_value);
      break;
    case PARILU_INPUT_TYPE:
      default_value.i = PARILU_TYPE;
      parilu_parse_number(&parilu->type, optarg, 0, default_value);
      break;
    case PARILU_INPUT_PIVOT:
      default_value.i = PARILU_PIVOT;
      parilu_parse_number(&parilu->pivot, optarg, 0, default_value);
      break;
    case PARILU_INPUT_TOL:
      default_value.d = PARILU_TOL;
      parilu_parse_number(&parilu->tol, optarg, 1, default_value);
      break;
    case PARILU_INPUT_NNZ_PER_ROW:
      default_value.i = PARILU_NNZ_PER_ROW;
      parilu_parse_number(&parilu->nnz_per_row, optarg, 0, default_value);
      break;
    case PARILU_INPUT_NULL_SPACE:
      parilu->null_space = atoi(optarg);
      break;
    case PARILU_INPUT_FILE:
      parilu->file = strndup(optarg, BUFSIZ);
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
struct parilu_opts_t *parilu_parse_opts(int *argc, char **argv[]) {
  struct parilu_opts_t *opts = parilu_calloc(struct parilu_opts_t, 1);

  parilu_parse_opts_aux(opts, argc, argv);

  return opts;
}

void parilu_finalize_opts(struct parilu_opts_t **opts) {
  parilu_free(&(*opts)->file);
  parilu_free(opts);
}
