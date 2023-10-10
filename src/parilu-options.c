#include <errno.h>
#include <getopt.h>

#include "parilu-impl.h"

/**
 * @ingroup parilu_user_api_functions
 *
 * @brief Get the default options for the parilu library. This can be used to
 * setup a system with parilu_setup(). User is responsible for calling
 * parilu_finalize_opts() to free the memory allocated by this function.
 *
 * @return ::parilu_opts*
 */
parilu_opts *parilu_default_opts() {
  parilu_opts *opts = parilu_calloc(parilu_opts, 1);
  return opts;
}

int parilu_set_verbose(parilu_opts *opts, unsigned verbose) {
  opts->verbose = verbose;
  return 0;
}

int parilu_set_matrix(parilu_opts *opts, const char *file) {
  opts->file = strndup(file, BUFSIZ);
  return 0;
}

int paril_set_null_space(parilu_opts *opts, unsigned null_space) {
  opts->null_space = null_space;
  return 0;
}

void parilu_finalize_opts(parilu_opts **opts) {
  parilu_free(&(*opts)->file);
  parilu_free(opts);
}
