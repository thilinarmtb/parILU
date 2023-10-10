#include <errno.h>
#include <getopt.h>

#include "parilu-impl.h"

/**
 * @ingroup parilu_user_api_functions
 *
 * @brief Get the default options for the parilu_handle library. This can be
 * used to setup a system with parilu_setup(). User is responsible for calling
 * parilu_finalize_opts() to free the memory allocated by this function.
 *
 * @return ::parilu_options*
 */
parilu_options *parilu_default_opts() {
  parilu_options *opts = parilu_calloc(parilu_options, 1);

  opts->verbose = PARILU_VERBOSE;
  opts->type = PARILU_TYPE;
  opts->pivot = PARILU_PIVOT;
  opts->tol = PARILU_TOL;
  opts->nnz_per_row = PARILU_NNZ_PER_ROW;
  opts->file = NULL;

  return opts;
}

int parilu_set_verbose(parilu_options *const opts, const unsigned verbose) {
  opts->verbose = verbose;
  return 0;
}

int parilu_get_verbose(const parilu_options *const opts, unsigned *verbose) {
  *verbose = opts->verbose;
  return 0;
}

int parilu_set_type(parilu_options *const opts, const unsigned type) {
  opts->type = type;
  return 0;
}

int parilu_get_type(const parilu_options *const opts, unsigned *type) {
  *type = opts->type;
  return 0;
}

int parilu_set_null_space(parilu_options *const opts,
                          const unsigned null_space) {
  opts->null_space = null_space;
  return 0;
}

int parilu_get_null_space(const parilu_options *const opts,
                          unsigned *null_space) {
  *null_space = opts->null_space;
  return 0;
}

int parilu_set_tol(parilu_options *const opts, const double tol) {
  opts->tol = tol;
  return 0;
}

int parilu_get_tol(const parilu_options *const opts, double *tol) {
  *tol = opts->tol;
  return 0;
}

int parilu_set_nnz_per_row(parilu_options *const opts,
                           const unsigned nnz_per_row) {
  opts->nnz_per_row = nnz_per_row;
  return 0;
}

int parilu_get_nnz_per_row(const parilu_options *const opts,
                           unsigned *nnz_per_row) {
  *nnz_per_row = opts->nnz_per_row;
  return 0;
}

int parilu_set_matrix(parilu_options *const opts, const char *const file) {
  opts->file = strndup(file, BUFSIZ);
  return 0;
}

int parilu_get_matrix(const parilu_options *const opts, char **const file) {
  *file = strndup(opts->file, BUFSIZ);
  return 0;
}

void parilu_finalize_opts(parilu_options **opts) {
  parilu_free(&(*opts)->file);
  parilu_free(opts);
}
