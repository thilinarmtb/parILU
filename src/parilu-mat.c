#include "parilu-impl.h"

#define MAX_ARRAY_SIZE 1024

struct mij_t {
  ulong r, c;
  uint p, idx;
  scalar v;
};

struct parilu_mat_t *parilu_mat_setup(const uint n, const slong *const vtx,
                                      const uint nnz, const uint *const row,
                                      const uint *const col,
                                      const double *const val,
                                      const struct comm *const c,
                                      buffer *const bfr, const int verbose) {
  parilu_debug(c, verbose, PARILU_INFO, "parilu_mat_setup: n = %u, nnz = %u.",
               n, nnz);

  // Check if the global matrix is empty and return if it is.
  struct parilu_mat_t *M = tcalloc(struct parilu_mat_t, 1);
  {
    slong nnzg = nnz, wrk;
    comm_allreduce(c, gs_long, gs_add, &nnzg, 1, &wrk);
    if (nnzg == 0) {
      M->rn = M->cn = 0;
      M->off = M->idx = NULL;
      M->row = M->col = NULL;
      return M;
    }
  }

  // Add the non-zeros into an array.
  struct array mat1;
  array_init(struct mij_t, &mat1, nnz);
  {
    slong minv = LONG_MAX;
    for (uint i = 0; i < nnz; i++) {
      parilu_assert(row[i] < n && col[i] < n,
                    "row[i] >= n and/or col[i] >= n.");
      ulong r = vtx[row[i]], c = vtx[col[i]];
      if (r == 0 || c == 0)
        continue;
      if ((slong)r < minv)
        minv = r;
      struct mij_t m = {.r = r, .c = c, .v = val[i]};
      array_cat(struct mij_t, &mat1, &m, 1);
    }

    slong wrk;
    comm_allreduce(c, gs_long, gs_min, &minv, 1, &wrk);
    parilu_assert(minv == 1, "minv != 1.");
  }

  // Assemble the entries locally by summing up the entries with same row and
  // column.
  struct array mat2;
  array_init(struct mij_t, &mat2, mat1.n);
  {
    sarray_sort_2(struct mij_t, mat1.ptr, mat1.n, r, 1, c, 1, bfr);
    struct mij_t *pm1 = (struct mij_t *)mat1.ptr;
    uint i = 0, j;
    while (i < mat1.n) {
      j = i + 1;
      while (j < mat1.n && (pm1[i].r == pm1[j].r) && (pm1[i].c == pm1[j].c))
        pm1[i].v += pm1[j].v, j++;
      pm1[i].p = pm1[i].r % c->np;
      array_cat(struct mij_t, &mat2, &pm1[i], 1);
      i = j;
    }
  }
  array_free(&mat1);

  // Assemble the matrix globally now to finish the assembly.
  uint rn = 0;
  {
    struct crystal cr;
    crystal_init(&cr, c);
    sarray_transfer(struct mij_t, &mat2, p, 1, &cr);
    crystal_free(&cr);

    array_init(struct mij_t, &mat1, mat2.n);
    if (mat2.n > 0) {
      sarray_sort_2(struct mij_t, mat2.ptr, mat2.n, r, 1, c, 1, bfr);
      struct mij_t *pm2 = (struct mij_t *)mat2.ptr;
      uint i = 0, j;
      while (i < mat2.n) {
        j = i + 1;
        while (j < mat2.n && (pm2[i].r == pm2[j].r) && (pm2[i].c == pm2[j].c))
          pm2[i].v += pm2[j].v, j++;
        array_cat(struct mij_t, &mat1, &pm2[i], 1);
        i = j, rn++;
      }
    }
  }
  array_free(&mat2);

  // Find the column ids in each processor.
  uint cn = 0, cmax = MAX_ARRAY_SIZE;
  ulong *cols = tcalloc(ulong, cmax);
  {
    sarray_sort(struct mij_t, mat1.ptr, mat1.n, c, 1, bfr);
    struct mij_t *pm1 = (struct mij_t *)mat1.ptr;
    uint i = 0, j;
    while (i < mat1.n) {
      pm1[i].idx = cn;
      j = i + 1;
      while (j < mat1.n && pm1[i].c == pm1[j].c)
        pm1[j].idx = cn, j++;

      if (cn == cmax) {
        cmax += cmax / 2 + 1;
        cols = trealloc(ulong, cols, cmax);
      }
      cols[cn] = pm1[i].c, cn++, i = j;
    }
  }

  // Setup the CSR matrix.
  {
    sarray_sort_2(struct mij_t, mat1.ptr, mat1.n, r, 1, c, 1, bfr);

    struct mij_t *pm1 = (struct mij_t *)mat1.ptr;
    M->rn = rn;
    M->off = tcalloc(uint, rn + 1);
    M->row = tcalloc(ulong, rn);
    M->idx = tcalloc(uint, mat1.n);
    M->val = tcalloc(scalar, mat1.n);

    uint i = 0, j, r = 0, nnz = 0;
    while (i < mat1.n) {
      M->idx[i] = pm1[i].idx, M->val[i] = pm1[i].v, nnz++;
      j = i + 1;
      while (j < mat1.n && pm1[i].r == pm1[j].r)
        M->idx[j] = pm1[j].idx, M->val[j] = pm1[j].v, j++, nnz++;
      M->row[r] = pm1[i].r, r++;
      M->off[r] = nnz;
    }
    // Check invariant: mat1.n == nnz
    parilu_assert(nnz == mat1.n, "nnz != mat1.n.");

    M->cn = cn;
    M->col = tcalloc(ulong, cn);
    for (uint i = 0; i < cn; i++)
      M->col[i] = cols[i];
  }
  parilu_free(&cols);

  return M;
}

struct parilu_mat_t *
parilu_mat_laplacian_setup(const struct parilu_mat_t *const M) {
  parilu_assert(M != NULL, "M == NULL.");
  parilu_assert(M->off != NULL, "M->off == NULL.");
  parilu_assert(M->idx != NULL, "M->idx == NULL.");
  parilu_assert(M->col != NULL, "M->col == NULL.");
  parilu_assert(M->row != NULL, "M->row == NULL.");

  struct parilu_mat_t *L = tcalloc(struct parilu_mat_t, 1);
  const uint rn = L->rn = M->rn;
  L->off = tcalloc(uint, rn + 1);
  memcpy(L->off, M->off, sizeof(uint) * (rn + 1));

  L->row = tcalloc(ulong, rn);
  memcpy(L->row, M->row, sizeof(ulong) * rn);

  const uint nnz = M->off[rn];
  L->idx = tcalloc(uint, nnz);
  memcpy(L->idx, M->idx, nnz);

  const uint cn = L->cn = M->cn;
  L->col = tcalloc(ulong, cn);
  memcpy(L->col, M->col, sizeof(ulong) * cn);

  for (uint i = 0; i < rn; i++) {
    sint didx = -1;
    for (uint j = M->idx[i], je = M->idx[i + 1]; j < je; j++) {
      if (M->row[i] == M->col[M->idx[j]])
        didx = j;
      L->val[j] = -1;
    }
    parilu_assert(didx != -1, "No diagonal entry found!");
    L->val[didx] = M->idx[i + 1] - M->idx[i] - 1;
  }

  return L;
}

void parilu_mat_free(struct parilu_mat_t **M_) {
  if (M_ && *M_) {
    struct parilu_mat_t *M = *M_;
    parilu_free(&M->off);
    parilu_free(&M->idx);
    parilu_free(&M->row);
    parilu_free(&M->col);
    parilu_free(&M->val);
  }
  parilu_free(M_);
}

struct parilu_mat_op_t *parilu_mat_op_setup(const struct parilu_mat_t *M,
                                            const struct comm *const c) {
  // Setup the gs handle for the communication required for the matrix-vector
  // product.
  const uint rn = M->rn, nnz = M->off[rn];

  // Initialize the buffer.
  struct parilu_mat_op_t *op = tcalloc(struct parilu_mat_op_t, 1);
  {
    buffer_init(&op->bfr, MAX_ARRAY_SIZE);
    buffer_reserve(&op->bfr, (nnz + rn) * sizeof(slong));
  }

  // Setup the gs handle.
  {
    slong *const ids = (slong *)op->bfr.ptr;
    for (uint i = 0; i < nnz; i++)
      ids[i] = -M->col[M->idx[i]];
    for (uint i = 0; i < rn; i++)
      ids[nnz + i] = M->row[i];
    op->gsh = gs_setup(ids, nnz + rn, c, 0, gs_pairwise, 0);
  }

  // Allocate work array set the pointer to the original matrix.
  {
    op->wrk = tcalloc(scalar, nnz + rn);
    op->M = M;
  }

  return op;
}

void parilu_mat_op_apply(scalar *const y, struct parilu_mat_op_t *const op,
                         const scalar *const x) {
  // Bring the entries necessry to do the mat-vec. Only rn entries are owned by
  // this process. Rest has to be brought in.
  const struct parilu_mat_t *const M = op->M;
  const uint rn = M->rn, nnz = M->off[rn];
  scalar *const wrk = op->wrk;
  {
    for (uint i = 0; i < rn; i++)
      wrk[nnz + i] = x[i];
    gs(wrk, gs_double, gs_add, 0, op->gsh, &op->bfr);
  }

  // Now perform the local mat-vec.
  {
    for (uint i = 0; i < rn; i++) {
      scalar sum = 0;
      for (uint j = M->idx[i], je = M->idx[i + 1]; j < je; j++)
        sum += M->val[j] * wrk[j];
      y[i] = sum;
    }
  }
}

void parilu_mat_op_free(struct parilu_mat_op_t **op_) {
  if (op_ && *op_) {
    struct parilu_mat_op_t *op = *op_;
    buffer_free(&op->bfr);
    gs_free(op->gsh);
    parilu_free(&op->wrk);
  }
  parilu_free(op_);
}

#undef MAX_ARRAY_SIZE
