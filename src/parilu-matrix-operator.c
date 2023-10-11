#include "parilu-impl.h"

parilu_matrix_operator *parilu_matrix_operator_setup(const parilu_matrix *M,
                                                     const MPI_Comm comm) {
  // Setup the gs handle for the communication required for the matrix-vector
  // product.
  const uint rn = M->rn, nnz = M->off[rn];

  // Initialize the buffer.
  parilu_matrix_operator *op = parilu_calloc(parilu_matrix_operator, 1);
  {
    buffer_init(&op->bfr, 1024);
    const size_t unit_size =
        sizeof(slong) > sizeof(scalar) ? sizeof(slong) : sizeof(scalar);
    buffer_reserve(&op->bfr, (nnz + rn) * unit_size);
  }

  // Setup the gs handle.
  {
    slong *const ids = (slong *const)op->bfr.ptr;
    for (uint i = 0; i < nnz; i++)
      ids[i] = -M->col[M->idx[i]];
    for (uint i = 0; i < rn; i++)
      ids[nnz + i] = M->row[i];

    struct comm c;
    comm_init(&c, comm);
    op->gsh = gs_setup(ids, nnz + rn, &c, 0, gs_pairwise, 0);
    comm_free(&c);
  }

  // Allocate work array set the pointer to the original matrix.
  {
    op->wrk = parilu_calloc(scalar, nnz + rn);
    op->M = M;
  }

  return op;
}

void parilu_matrix_operator_apply(scalar *const y,
                                  parilu_matrix_operator *const op,
                                  const scalar *const x) {
  // Bring the entries necessry to do the mat-vec. Only rn entries are owned by
  // this process. Rest has to be brought in.
  const parilu_matrix *const M = op->M;
  const uint rn = M->rn;
  scalar *const wrk = (scalar *const)op->wrk;
  {
    const uint nnz = M->off[rn];
    for (uint i = 0; i < nnz; i++)
      wrk[i] = 0;
    for (uint i = 0; i < rn; i++)
      wrk[nnz + i] = x[i];
    gs(wrk, gs_double, gs_add, 0, op->gsh, &op->bfr);
  }

  // Now perform the local mat-vec.
  {
    for (uint i = 0; i < rn; i++) {
      scalar sum = 0;
      for (uint j = M->off[i], je = M->off[i + 1]; j < je; j++)
        sum += M->val[j] * wrk[j];
      y[i] = sum;
    }
  }
}

void parilu_matrix_operator_free(parilu_matrix_operator **op_) {
  if (op_ && *op_) {
    parilu_matrix_operator *op = *op_;
    buffer_free(&op->bfr);
    gs_free(op->gsh);
    parilu_free(&op->wrk);
  }
  parilu_free(op_);
}
