#include "parilu-impl.h"

#define MAX_ARRAY_SIZE 1024

struct mij_t {
  ulong r, c;
  uint p, idx;
  scalar v;
};

parilu_matrix *parilu_matrix_setup(const uint32_t nnz,
                                   const uint64_t *const row,
                                   const uint64_t *const col,
                                   const double *const val,
                                   const MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  // Check if the global matrix is empty and return if it is.
  parilu_matrix *M = parilu_calloc(parilu_matrix, 1);
  {
    slong nnzg = nnz, wrk;
    comm_allreduce(&c, gs_long, gs_add, &nnzg, 1, &wrk);
    parilu_log(&c, PARILU_INFO, "parilu_matrix_setup: nnzg = %ld.", nnzg);
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
    comm_allreduce(&c, gs_long, gs_min, &minv, 1, &wrk);
    parilu_log(&c, PARILU_INFO, "parilu_matrix_setup: minv = %ld.", minv);
    parilu_assert(minv == 1, "minv != 1");
  }

  // Assemble the entries locally by summing up the entries with same row and
  // column.
  struct array mat2;
  buffer bfr;
  {
    buffer_init(&bfr, 1024);
    sarray_sort_2(struct mij_t, mat1.ptr, mat1.n, r, 1, c, 1, &bfr);
    struct mij_t *pm1 = (struct mij_t *)mat1.ptr;

    array_init(struct mij_t, &mat2, mat1.n);
    uint i = 0, j;
    while (i < mat1.n) {
      j = i + 1;
      while (j < mat1.n && (pm1[i].r == pm1[j].r) && (pm1[i].c == pm1[j].c))
        pm1[i].v += pm1[j].v, j++;
      pm1[i].p = pm1[i].r % c.np;
      array_cat(struct mij_t, &mat2, &pm1[i], 1);
      i = j;
    }
    array_free(&mat1);
  }

  // Assemble the matrix globally now to finish the assembly.
  {
    struct crystal cr;
    crystal_init(&cr, &c);
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
  {
    sarray_sort(struct mij_t, mat1.ptr, mat1.n, c, 1, &bfr);
    struct mij_t *pm1 = (struct mij_t *)mat1.ptr;

    uint cn = 0, cmax = MAX_ARRAY_SIZE;
    M->col = parilu_calloc(ulong, cmax);

    uint i = 0;
    while (i < mat1.n) {
      pm1[i].idx = cn;
      uint j = i + 1;
      while (j < mat1.n && pm1[i].c == pm1[j].c)
        pm1[j].idx = cn, j++;

      if (cn == cmax) {
        cmax += cmax / 2 + 1;
        M->col = trealloc(ulong, M->col, cmax);
      }
      M->col[cn] = pm1[i].c, cn++, i = j;
    }

    M->cn = cn;
  }

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
    M->off[0] = 0;
    while (i < mat1.n) {
      uint j = i;
      while (j < mat1.n && pm1[i].r == pm1[j].r)
        M->idx[j] = pm1[j].idx, M->val[j] = pm1[j].v, j++, nnz++;
      M->row[r] = pm1[i].r, r++;
      M->off[r] = nnz;
      i = j;
    }

    // Check invariant: r == rn
    parilu_assert(r == rn, "r != rn.");
    // Check invariant: mat1.n == nnz
    parilu_assert(nnz == mat1.n, "nnz != mat1.n.");

    buffer_free(&bfr), array_free(&mat1);
  }
  parilu_log(&c, PARILU_INFO, "parilu_matrix_setup: done.");

  comm_free(&c);

  return M;
}

parilu_matrix *parilu_matrix_laplacian_setup(const parilu_matrix *const M) {
  parilu_assert(M != NULL, "M == NULL.");
  parilu_assert(M->rn > 0, "M->rn == 0");

  parilu_matrix *L = parilu_calloc(parilu_matrix, 1);

  const uint rn = L->rn = M->rn;
  L->off = parilu_calloc(uint, rn + 1);
  memcpy(L->off, M->off, sizeof(uint) * (rn + 1));

  L->row = parilu_calloc(ulong, rn);
  memcpy(L->row, M->row, sizeof(ulong) * rn);

  const uint nnz = M->off[rn];
  L->idx = parilu_calloc(uint, nnz);
  memcpy(L->idx, M->idx, sizeof(ulong) * nnz);

  const uint cn = L->cn = M->cn;
  L->col = parilu_calloc(ulong, cn);
  memcpy(L->col, M->col, sizeof(ulong) * cn);

  L->val = parilu_calloc(scalar, nnz);
  for (uint i = 0; i < rn; i++) {
    sint didx = -1;
    for (uint j = M->off[i], je = M->off[i + 1]; j < je; j++) {
      if (M->row[i] == M->col[M->idx[j]])
        didx = j;
      L->val[j] = -1;
    }
    parilu_assert(didx >= 0, "No diagonal entry found!");
    L->val[didx] = (scalar)(M->off[i + 1] - M->off[i]) - 1;
  }

  return L;
}

void parilu_matrix_dump(const char *const file, const parilu_matrix *const M,
                        const MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  parilu_log(&c, PARILU_INFO, "parilu_matrix_dump: file = %s.", file);

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
    crystal_init(&cr, &c);
    sarray_transfer(struct data_t, &arr, p, 0, &cr);
    crystal_free(&cr);

    buffer bfr;
    buffer_init(&bfr, 1024);
    sarray_sort_2(struct data_t, arr.ptr, arr.n, r, 1, c, 1, &bfr);
    buffer_free(&bfr);
  }

  // Write the data to file at rank 0.
  if (c.id == 0) {
    FILE *fp = fopen(file, "w+");
    if (!fp) {
      fprintf(stderr, "Failed to open file %s.\n", file);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const struct data_t *const ptr = (const struct data_t *)arr.ptr;
    for (uint i = 0; i < arr.n; i++)
      fprintf(fp, "%llu %llu %e\n", ptr[i].r, ptr[i].c, ptr[i].v);

    fclose(fp);
  }

  array_free(&arr);
  comm_free(&c);
}

void parilu_matrix_free(parilu_matrix **M_) {
  if (M_ && *M_) {
    parilu_matrix *M = *M_;
    parilu_free(&M->off);
    parilu_free(&M->idx);
    parilu_free(&M->row);
    parilu_free(&M->col);
    parilu_free(&M->val);
  }
  parilu_free(M_);
}

#undef MAX_ARRAY_SIZE
