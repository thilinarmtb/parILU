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

  opts->verbose = PARILU_VERBOSE;
  opts->type = PARILU_TYPE;
  opts->pivot = PARILU_PIVOT;
  opts->tol = PARILU_TOL;
  opts->nnz_per_row = PARILU_NNZ_PER_ROW;
  opts->file = NULL;

  return opts;
}

int parilu_set_verbose(parilu_opts *const opts, const unsigned verbose) {
  opts->verbose = verbose;
  return 0;
}

int parilu_get_verbose(const parilu_opts *const opts, unsigned *verbose) {
  *verbose = opts->verbose;
  return 0;
}

int parilu_set_type(parilu_opts *const opts, const unsigned type) {
  opts->type = type;
  return 0;
}

int parilu_get_type(const parilu_opts *const opts, unsigned *type) {
  *type = opts->type;
  return 0;
}

int parilu_set_null_space(parilu_opts *const opts, const unsigned null_space) {
  opts->null_space = null_space;
  return 0;
}

int parilu_get_null_space(const parilu_opts *const opts, unsigned *null_space) {
  *null_space = opts->null_space;
  return 0;
}

int parilu_set_tol(parilu_opts *const opts, const double tol) {
  opts->tol = tol;
  return 0;
}

int parilu_get_tol(const parilu_opts *const opts, double *tol) {
  *tol = opts->tol;
  return 0;
}

int parilu_set_nnz_per_row(parilu_opts *const opts,
                           const unsigned nnz_per_row) {
  opts->nnz_per_row = nnz_per_row;
  return 0;
}

int parilu_get_nnz_per_row(const parilu_opts *const opts,
                           unsigned *nnz_per_row) {
  *nnz_per_row = opts->nnz_per_row;
  return 0;
}

int parilu_set_matrix(parilu_opts *const opts, const char *const file) {
  opts->file = strndup(file, BUFSIZ);
  return 0;
}

int parilu_get_matrix(const parilu_opts *const opts, char **const file) {
  *file = strndup(opts->file, BUFSIZ);
  return 0;
}

void parilu_finalize_opts(parilu_opts **opts) {
  parilu_free(&(*opts)->file);
  parilu_free(opts);
}
