#include <errno.h>
#include <getopt.h>

#include "parilu-impl.h"

/**
 * @ingroup parilu_user_api_functions
 *
 * @brief Get the default options for the parilu_handle library. This can be
 * used to setup a system with parilu_setup(). User is responsible for calling
 * parilu_finalize_options() to free the memory allocated by this function.
 *
 * @return ::parilu_options*
 */
parilu_options *parilu_default_options() {
  parilu_options *options = parilu_calloc(parilu_options, 1);

  options->type = PARILU_TYPE;
  options->pivot = PARILU_PIVOT;
  options->tol = PARILU_TOL;
  options->nnz_per_row = PARILU_NNZ_PER_ROW;

  return options;
}

int parilu_set_type(parilu_options *const options, const unsigned type) {
  options->type = type;
  return 0;
}

int parilu_get_type(const parilu_options *const options, unsigned *type) {
  *type = options->type;
  return 0;
}

int parilu_set_null_space(parilu_options *const options,
                          const unsigned null_space) {
  options->null_space = null_space;
  return 0;
}

int parilu_get_null_space(const parilu_options *const options,
                          unsigned *null_space) {
  *null_space = options->null_space;
  return 0;
}

int parilu_set_tol(parilu_options *const options, const double tol) {
  options->tol = tol;
  return 0;
}

int parilu_get_tol(const parilu_options *const options, double *tol) {
  *tol = options->tol;
  return 0;
}

int parilu_set_nnz_per_row(parilu_options *const options,
                           const unsigned nnz_per_row) {
  options->nnz_per_row = nnz_per_row;
  return 0;
}

int parilu_get_nnz_per_row(const parilu_options *const options,
                           unsigned *nnz_per_row) {
  *nnz_per_row = options->nnz_per_row;
  return 0;
}

void parilu_finalize_options(parilu_options **options) { parilu_free(options); }
