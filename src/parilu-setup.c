#include "parilu-impl.h"

/**
 * @ingroup parilu_user_api_functions
 * @brief Setup parilu. Returns a pointer to a newly allocated struct parilu_t.
 *
 * @param n Number of local dofs in vertex array.
 * @param vertex Array of local dofs with global numbering.
 * @param nnz Number of nonzeros in the matrix.
 * @param row Row indices of the matrix which points to a global dof in \p
 * vertex array.
 * @param col Column indices of the matrix which points to a global dof in \p
 * vertex array.
 * @param val Values of the matrix.
 * @param options Pointer to the struct parilu_opts_t which contains the
 * options.
 * @param comm MPI communicator.
 * @param bfr Pointer to the buffer struct used for work arrays.
 */
struct parilu_t *parilu_setup(const uint n, const slong *const vertex,
                              const uint nnz, const uint *const row,
                              const uint *const col, const double *const val,
                              const struct parilu_opts_t *const options,
                              const MPI_Comm comm, buffer *const bfr) {
  const int verbose = options->verbose;
  parilu_log_init(verbose);

  // Create a gslib comm out of MPI_Comm
  struct comm c;
  comm_init(&c, comm);

  // Initialize ILU struct.
  parilu_log(&c, PARILU_INFO, "parilu_setup: Initialize ILU options.");
  struct parilu_t *ilu = tcalloc(struct parilu_t, 1);
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
  parilu_log(&c, PARILU_INFO, "parilu_setup: Setup the matrix.");
  struct parilu_mat_t *M =
      parilu_mat_setup(n, vertex, nnz, row, col, val, &c, bfr, verbose - 1);

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
