#include "parilu-impl.h"

/**
 * @ingroup parilu_user_api_functions
 * @brief Setup parilu. Returns a pointer to a newly allocated struct parilu_t.
 *
 * @param nnz Number of nonzeros in the matrix.
 * @param row Global row index or number.
 * @param col Global column index or number.
 * @param val Values of the matrix. val[i] is the value of the matrix entry
 * (row[i], col[i]).
 * @param options Pointer to the struct parilu_opts_t which contains the
 * options.
 * @param comm MPI communicator.
 */
struct parilu_t *parilu_setup(const uint32_t nnz, const uint64_t *const row,
                              const uint64_t *const col,
                              const double *const val,
                              const struct parilu_opts_t *const options,
                              const MPI_Comm comm) {
  // Create a gslib comm out of MPI_Comm.
  struct comm c;
  comm_init(&c, comm);

  const int verbose = options->verbose;
  parilu_set_log_level(verbose);

  parilu_log(&c, PARILU_INFO, "parilu_setup: verbose = %d", options->verbose);
  parilu_log(&c, PARILU_INFO, "parilu_setup: type = %d", options->type);
  parilu_log(&c, PARILU_INFO, "parilu_setup: pivot = %d", options->pivot);
  parilu_log(&c, PARILU_INFO, "parilu_setup: tol = %e", options->tol);
  parilu_log(&c, PARILU_INFO, "parilu_setup: nnz_per_row = %d",
             options->nnz_per_row);
  parilu_log(&c, PARILU_INFO, "parilu_setup: null_space = %d",
             options->null_space);
  parilu_log(&c, PARILU_INFO, "parilu_setup: file = %s", options->file);

  // Initialize ILU struct.
  parilu_log(&c, PARILU_INFO, "parilu_setup: Initialize ILU options.");
  struct parilu_t *ilu = parilu_calloc(struct parilu_t, 1);
  {
    ilu->pivot = options->pivot;
    ilu->null_space = options->null_space;
    ilu->tol = options->tol;
    ilu->nnz_per_row = options->nnz_per_row;
    // Things to be set during factorization.
    ilu->nlvls = 0;
    ilu->lvl_off = NULL;
    ilu->perm = NULL;
  }

  // Setup CSR mat for ILU system.
  parilu_log(&c, PARILU_INFO, "parilu_setup: Setup the CSR matrix.");
  struct parilu_mat_t *M =
      parilu_mat_setup(nnz, row, col, val, &c, verbose - 1);

  parilu_mat_dump("system.txt", M, &c);

  // Create the Laplacian matrix of the system.
  parilu_log(&c, PARILU_INFO, "parilu_setup: Setup the Laplacian matrix.");
  struct parilu_mat_t *L = parilu_mat_laplacian_setup(M);

  // Parition the matrix with parRSB.
  parilu_log(&c, PARILU_INFO, "parilu_setup: Partition the matrix.");

  parilu_mat_free(&M), parilu_mat_free(&L);
  comm_free(&c);

  parilu_log(&c, PARILU_INFO, "parilu_setup: done.");

  return ilu;
}
