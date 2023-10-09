#include "parilu-impl.h"

#define MAX_ARRAY_SIZE 1024

struct mij_t {
  ulong r, c;
  uint p, idx;
  scalar v;
};

struct parilu_mat_t *
parilu_mat_setup(const uint32_t nnz, const uint64_t *const row,
                 const uint64_t *const col, const double *const val,
                 const struct comm *const c, const int verbose) {
  // Check if the global matrix is empty and return if it is.
  struct parilu_mat_t *M = parilu_calloc(struct parilu_mat_t, 1);
  {
    slong nnzg = nnz, wrk;
    comm_allreduce(c, gs_long, gs_add, &nnzg, 1, &wrk);
    parilu_log(c, PARILU_INFO, "parilu_mat_setup: nnzg = %ld.", nnzg);
    if (nnzg == 0) {
      M->rn = M->cn = 0;
      M->off = M->idx = NULL;
      M->row = M->col = NULL;
      return M;
    }
  }

  // Add the non-zeros into an array.
  struct array mat1;
  {
    array_init(struct mij_t, &mat1, nnz);
    slong minv = LONG_MAX;
    for (uint32_t i = 0; i < nnz; i++) {
      ulong r = row[i], c = col[i];
      if (r == 0 || c == 0)
        continue;
      if ((slong)r < minv)
        minv = r;
      struct mij_t m = {.r = r, .c = c, .v = val[i]};
      array_cat(struct mij_t, &mat1, &m, 1);
    }

    slong wrk;
    comm_allreduce(c, gs_long, gs_min, &minv, 1, &wrk);
    parilu_log(c, PARILU_INFO, "parilu_mat_setup: minv = %ld.", minv);
    parilu_assert(minv == 1, "minv != 1");
  }

  buffer bfr;
  buffer_init(&bfr, 1024);

  // Assemble the entries locally by summing up the entries with same row and
  // column.
  struct array mat2;
  {
    sarray_sort_2(struct mij_t, mat1.ptr, mat1.n, r, 1, c, 1, &bfr);
    struct mij_t *pm1 = (struct mij_t *)mat1.ptr;

    array_init(struct mij_t, &mat2, mat1.n);
    uint i = 0, j;
    while (i < mat1.n) {
      j = i + 1;
      while (j < mat1.n && (pm1[i].r == pm1[j].r) && (pm1[i].c == pm1[j].c))
        pm1[i].v += pm1[j].v, j++;
      pm1[i].p = pm1[i].r % c->np;
      array_cat(struct mij_t, &mat2, &pm1[i], 1);
      i = j;
    }
    array_free(&mat1);
  }

  // Assemble the matrix globally now to finish the assembly.
  {
    struct crystal cr;
    crystal_init(&cr, c);
    sarray_transfer(struct mij_t, &mat2, p, 1, &cr);
    crystal_free(&cr);

    sarray_sort_2(struct mij_t, mat2.ptr, mat2.n, r, 1, c, 1, &bfr);
    struct mij_t *pm2 = (struct mij_t *)mat2.ptr;

    array_init(struct mij_t, &mat1, mat2.n);
    uint i = 0, j;
    while (i < mat2.n) {
      j = i + 1;
      while (j < mat2.n && (pm2[i].r == pm2[j].r) && (pm2[i].c == pm2[j].c))
        pm2[i].v += pm2[j].v, j++;
      array_cat(struct mij_t, &mat1, &pm2[i], 1);
      i = j;
    }
    array_free(&mat2);
  }
  parilu_log(c, PARILU_INFO, "parilu_mat_setup: nnz = %u.", mat1.n);

  // Find the number of rows in each processor.
  uint rn = 0;
  {
    const struct mij_t *const pm1 = (const struct mij_t *const)mat1.ptr;
    uint i = 0;
    while (i < mat1.n) {
      uint j = i + 1;
      while (j < mat1.n && pm1[i].r == pm1[j].r)
        j++;
      rn++, i = j;
    }
  }

  // Find the column ids in each processor.
  uint cn = 0, cmax = MAX_ARRAY_SIZE;
  ulong *cols = parilu_calloc(ulong, cmax);
  {
    sarray_sort(struct mij_t, mat1.ptr, mat1.n, c, 1, &bfr);
    struct mij_t *pm1 = (struct mij_t *)mat1.ptr;

    uint i = 0;
    while (i < mat1.n) {
      pm1[i].idx = cn;
      uint j = i + 1;
      while (j < mat1.n && pm1[i].c == pm1[j].c)
        pm1[j].idx = cn, j++;

      if (cn == cmax) {
        cmax += cmax / 2 + 1;
        cols = trealloc(ulong, cols, cmax);
      }
      cols[cn] = pm1[i].c, cn++, i = j;
    }
  }
  parilu_log(c, PARILU_INFO, "parilu_mat_setup: cn = %u.", cn);

  // Setup the CSR matrix.
  {
    sarray_sort_2(struct mij_t, mat1.ptr, mat1.n, r, 1, c, 1, &bfr);
    struct mij_t *pm1 = (struct mij_t *)mat1.ptr;

    M->rn = rn;
    M->off = parilu_calloc(uint, rn + 1);
    M->row = parilu_calloc(ulong, rn);
    M->idx = parilu_calloc(uint, mat1.n);
    M->val = parilu_calloc(scalar, mat1.n);

    uint i = 0, r = 0, nnz = 0;
    while (i < mat1.n) {
      uint j = i;
      while (j < mat1.n && pm1[i].r == pm1[j].r)
        M->idx[j] = pm1[j].idx, M->val[j] = pm1[j].v, j++, nnz++;
      M->row[r] = pm1[i].r, r++;
      M->off[r] = nnz;
      i = j;
    }

    parilu_log(c, PARILU_INFO, "parilu_mat_setup: r = %u.", r);
    // Check invariant: r == rn
    parilu_assert(r == rn, "r != rn.");
    // Check invariant: mat1.n == nnz
    parilu_assert(nnz == mat1.n, "nnz != mat1.n.");
    array_free(&mat1);

    M->cn = cn;
    M->col = parilu_calloc(ulong, cn);
    for (uint i = 0; i < cn; i++)
      M->col[i] = cols[i];

    parilu_free(&cols);
  }
  parilu_log(c, PARILU_INFO, "parilu_mat_setup: nnz = %u.", mat1.n);

  buffer_free(&bfr);

  return M;
}

struct parilu_mat_t *
parilu_mat_laplacian_setup(const struct parilu_mat_t *const M) {
  parilu_assert(M != NULL, "M == NULL.");
  parilu_assert(M->off != NULL, "M->off == NULL.");
  parilu_assert(M->idx != NULL, "M->idx == NULL.");
  parilu_assert(M->col != NULL, "M->col == NULL.");
  parilu_assert(M->row != NULL, "M->row == NULL.");

  struct parilu_mat_t *L = parilu_calloc(struct parilu_mat_t, 1);
  const uint rn = L->rn = M->rn;
  L->off = parilu_calloc(uint, rn + 1);
  memcpy(L->off, M->off, sizeof(uint) * (rn + 1));

  L->row = parilu_calloc(ulong, rn);
  memcpy(L->row, M->row, sizeof(ulong) * rn);

  const uint nnz = M->off[rn];
  L->idx = parilu_calloc(uint, nnz);
  memcpy(L->idx, M->idx, nnz);

  const uint cn = L->cn = M->cn;
  L->col = parilu_calloc(ulong, cn);
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

PARILU_INTERN void parilu_mat_dump(const char *const file,
                                   const struct parilu_mat_t *const M,
                                   const struct comm *const c) {
  struct data_t {
    ulong r, c;
    uint p;
    scalar v;
  };

  // Copy the matrix data to an array.
  struct array arr;
  {
    array_init(struct data_t, &arr, M->off[M->rn]);

    struct data_t dt = {.p = 0};

    for (uint i = 0; i < M->rn; i++) {
      for (uint j = M->off[i], je = M->off[i + 1]; j < je; j++) {
        dt.r = M->row[i], dt.c = M->col[M->idx[j]], dt.v = M->val[j];
        array_cat(struct data_t, &arr, &dt, 1);
      }
    }
  }

  // Send the data to rank 0.
  {
    struct crystal cr;
    crystal_init(&cr, c);
    sarray_transfer(struct data_t, &arr, p, 0, &cr);
    crystal_free(&cr);

    buffer bfr;
    buffer_init(&bfr, 1024);
    sarray_sort_2(struct data_t, arr.ptr, arr.n, r, 1, c, 1, &bfr);
    buffer_free(&bfr);
  }

  // Free data and exit if not rank 0.
  if (c->id) {
    array_free(&arr);
    return;
  }

  // Write the data to file at rank 0.
  {
    FILE *fp = fopen(file, "w");
    if (!fp)
      parilu_log(c, PARILU_ERROR, "Failed to open file %s.", file);

    const struct data_t *const ptr = (const struct data_t *)arr.ptr;
    for (uint i = 0; i < arr.n; i++)
      fprintf(fp, "%llu %llu %e\n", ptr[i].r, ptr[i].c, ptr[i].v);

    array_free(&arr);
    fclose(fp);
  }
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
  struct parilu_mat_op_t *op = parilu_calloc(struct parilu_mat_op_t, 1);
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
    op->wrk = parilu_calloc(scalar, nnz + rn);
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
