#include "parilu-impl.h"

/**
 * @ingroup parilu_user_api_functions
 * @brief Setup parilu. Returns a pointer to a newly allocated ::parilu_handle
 * handle.
 *
 * @param nnz Number of nonzeros in the matrix.
 * @param row Global row index or number.
 * @param col Global column index or number.
 * @param val Values of the matrix. val[i] is the value of the matrix entry
 * (row[i], col[i]).
 * @param options Pointer to the ::parilu_options which stores the options.
 * @param comm MPI communicator.
 */
parilu_handle *parilu_setup(const uint32_t nnz, const uint64_t *const row,
                            const uint64_t *const col, const double *const val,
                            const parilu_options *const options,
                            const MPI_Comm comm) {
  // Create a gslib comm out of MPI_Comm.
  struct comm c;
  comm_init(&c, comm);

  parilu_log(&c, PARILU_INFO, "parilu_setup: type = %d", options->type);
  parilu_log(&c, PARILU_INFO, "parilu_setup: pivot = %d", options->pivot);
  parilu_log(&c, PARILU_INFO, "parilu_setup: tol = %e", options->tol);
  parilu_log(&c, PARILU_INFO, "parilu_setup: nnz_per_row = %d",
             options->nnz_per_row);
  parilu_log(&c, PARILU_INFO, "parilu_setup: null_space = %d",
             options->null_space);

  // Initialize ILU struct.
  parilu_log(&c, PARILU_INFO, "parilu_setup: Initialize ILU options.");
  parilu_handle *ilu = parilu_calloc(parilu_handle, 1);
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
  parilu_matrix *M = parilu_matrix_setup(nnz, row, col, val, comm);

  // Create the Laplacian matrix of the system.
  parilu_log(&c, PARILU_INFO, "parilu_matrix_laplacian_setup: ...");
  parilu_matrix *L = parilu_matrix_laplacian_setup(M);

  // Parition the matrix with parRSB.
  parilu_partition(L, comm);
  parilu_matrix_free(&L);
  parilu_matrix_free(&M);

  parilu_log(&c, PARILU_INFO, "parilu_setup: done.");
  comm_free(&c);

  return ilu;
}
