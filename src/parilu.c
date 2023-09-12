#include <getopt.h>
#include <parilu-impl.h>

static void parilu_print_help(const char *name) {
  printf("Usage: %s [OPTIONS]\n", name);
  printf("Options:\n");
  printf("  --parilu-verbose=<verbose level>, Verbose level (0, 1, 2, ...).\n");
  printf("  --parilu-help, Prints this help message and exit.\n");
}

static void parilu_parse_opts(struct parilu_opts_t *parilu, int *argc,
                              char ***argv_) {
  parilu->verbose = 0;

  enum parilu_input_t { PARILU_INPUT_VERBOSE = 0, PARILU_INPUT_HELP = 99 };
  static struct option long_options[] = {
      {"parilu-verbose", optional_argument, 0, PARILU_INPUT_VERBOSE},
      {"parilu-help", no_argument, 0, PARILU_INPUT_HELP},
      {0, 0, 0, 0}};

  char *end = NULL;
  char **argv = *argv_;
  for (;;) {
    int opt = getopt_long(*argc, argv, "", long_options, NULL);
    if (opt == -1)
      break;

    switch (opt) {
    case PARILU_INPUT_VERBOSE:
      // NOLINTNEXTLINE
      parilu->verbose = strtol(optarg, &end, 10);
      if (!end || *end != '\0' || optarg == end)
        parilu_error("Invalid string for verbose level: %s\n", optarg);
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
 * @brief Initialize parilu. Returns a pointer to a newly allocated struct
 * parilu_opts_t.
 *
 * @param argc Pointer to the number of commnad line arguments.
 * @param argv Pointer to the array of command line arguments.
 *
 * @return struct parilu_opts_t*
 */
struct parilu_opts_t *parilu_init(int *argc, char **argv[]) {
  struct parilu_opts_t *parilu = parilu_calloc(struct parilu_opts_t, 1);

  parilu_parse_opts(parilu, argc, argv);

  return parilu;
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
