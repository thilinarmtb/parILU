#include "parilu-impl.h"

#define MAX_ARRAY_SIZE 1024

struct mij_t {
  ulong r, c;
  uint p, idx;
  scalar v;
};

struct csr_mat {
  uint n, cn, *off, *idx;
  ulong *row, *col;
  scalar *val;
};

static void parilu_mat_free(struct csr_mat *const M) {
  if (M) {
    if (M->off)
      free(M->off);
    if (M->idx)
      free(M->idx);
    if (M->row)
      free(M->row);
    if (M->col)
      free(M->col);
    if (M->val)
      free(M->val);
    free(M);
  }
}

static struct csr_mat *parilu_setup_mat(uint n, const slong *const vtx,
                                        const uint nnz, const uint *const row,
                                        const uint *const col,
                                        const double *const val,
                                        const struct comm *const c,
                                        buffer *bfr) {
  // Check if the global matrix is empty and return if it is.
  struct csr_mat *M = tcalloc(struct csr_mat, 1);
  {
    slong nnzg = nnz, wrk;
    comm_allreduce(c, gs_long, gs_add, &nnzg, 1, &wrk);
    if (nnzg == 0) {
      M->n = 0;
      M->off = M->idx = NULL;
      M->row = M->col = NULL;
      return M;
    }
  }

  // Add the non-zeros into an array.
  struct array mat;
  array_init(struct mij_t, &mat, nnz);
  {
    slong minv = LONG_MAX;
    for (uint i = 0; i < nnz; i++) {
      assert(row[i] < n && col[i] < n && "row[i] >= n and/or col[i] >= n");
      ulong r = vtx[row[i]], c = vtx[col[i]];
      if (r == 0 || c == 0)
        continue;
      if (r < (ulong)minv)
        minv = r;
      struct mij_t m = {.r = r, .c = c, .v = val[i]};
      array_cat(struct mij_t, &mat, &m, 1);
    }

    slong wrk;
    comm_allreduce(c, gs_long, gs_min, &minv, 1, &wrk);
    assert(minv == 1 && "minv != 1");
  }

  // Assemble the entries locally by summing up the entries with same row and
  // column.
  struct array umat;
  array_init(struct mij_t, &umat, mat.n);
  {
    sarray_sort_2(struct mij_t, mat.ptr, mat.n, r, 1, c, 1, bfr);
    struct mij_t *pm = (struct mij_t *)mat.ptr;
    uint i = 0, j;
    while (i < mat.n) {
      j = i + 1;
      while (j < mat.n && (pm[i].r == pm[j].r) && (pm[i].c == pm[j].c))
        pm[i].v += pm[j].v, j++;
      pm[i].p = pm[i].r % c->np;
      array_cat(struct mij_t, &umat, &pm[i], 1);
      i = j;
    }
    array_free(&mat);
  }

  // Assemble the matrix globally now to finish the global assembly.
  uint rn = 0;
  {
    struct crystal cr;
    crystal_init(&cr, c);
    sarray_transfer(struct mij_t, &umat, p, 1, &cr);
    crystal_free(&cr);

    array_init(struct mij_t, &mat, umat.n);
    if (umat.n > 0) {
      sarray_sort_2(struct mij_t, umat.ptr, umat.n, r, 1, c, 1, bfr);
      struct mij_t *pu = (struct mij_t *)umat.ptr;
      uint i = 0, j;
      while (i < umat.n) {
        j = i + 1;
        while (j < umat.n && (pu[i].r == pu[j].r) && (pu[i].c == pu[j].c))
          pu[i].v += pu[j].v, j++;
        array_cat(struct mij_t, &mat, &pu[i], 1);
        i = j, rn++;
      }
    }
    array_free(&umat);
  }

  // Find the column ids in each processor.
  uint cn = 0;
  ulong *cols = NULL;
  {
    uint cmax = MAX_ARRAY_SIZE;
    cols = tcalloc(ulong, cmax);

    sarray_sort(struct mij_t, mat.ptr, mat.n, c, 1, bfr);
    struct mij_t *pm = (struct mij_t *)mat.ptr;
    uint i = 0, j;
    while (i < mat.n) {
      pm[i].idx = cn;
      j = i + 1;
      while (j < mat.n && pm[i].c == pm[j].c)
        pm[j].idx = cn, j++;

      if (cn == cmax) {
        cmax += cmax / 2 + 1;
        cols = trealloc(ulong, cols, cmax);
      }
      cols[cn] = pm[i].c, cn++, i = j;
    }
  }

  // Setup the CSR matrix.
  {
    sarray_sort_2(struct mij_t, mat.ptr, mat.n, r, 1, c, 1, bfr);

    struct mij_t *pm = (struct mij_t *)mat.ptr;
    M->n = rn;
    M->off = tcalloc(uint, rn + 1);
    M->row = tcalloc(ulong, rn);
    M->idx = tcalloc(uint, mat.n);
    M->val = tcalloc(scalar, mat.n);

    uint i = 0, j, r = 0, nnz = 0;
    while (i < mat.n) {
      M->idx[i] = pm[i].idx, M->val[i] = pm[i].v, nnz++;
      j = i + 1;
      while (j < mat.n && pm[i].r == pm[j].r)
        M->idx[j] = pm[j].idx, M->val[j] = pm[j].v, j++, nnz++;
      M->row[r] = pm[i].r, r++;
      M->off[r] = nnz;
    }
    // Check invariant: mat.n == nnz
    assert(nnz == mat.n);

    M->cn = cn;
    M->col = tcalloc(ulong, cn);
    for (uint i = 0; i < cn; i++)
      M->col[i] = cols[i];
    free(cols);
  }

  return M;
}

static struct csr_mat *
parilu_setup_laplacian_mat(const struct csr_mat *const M) {
  struct csr_mat *L = tcalloc(struct csr_mat, 1);

  const uint n = L->n = M->n;
  L->off = tcalloc(uint, M->n + 1);
  memcpy(L->off, M->off, sizeof(uint) * (n + 1));

  L->row = tcalloc(ulong, n);
  memcpy(L->row, M->row, sizeof(ulong) * n);

  const uint nnz = M->off[n];
  L->idx = tcalloc(uint, nnz);
  memcpy(L->idx, M->idx, nnz);

  const uint cn = L->cn = M->cn;
  L->col = tcalloc(ulong, cn);
  memcpy(L->col, M->col, sizeof(ulong) * cn);

  for (uint i = 0; i < n; i++) {
    sint didx = -1;
    for (uint j = M->idx[i], je = M->idx[i + 1]; j < je; j++) {
      if (M->row[i] == M->col[M->idx[j]])
        didx = j;
      M->val[j] = -1;
    }
    assert(didx >= M->idx[i] && didx < M->idx[i + 1]);
    assert((M->idx[i + 1] - M->idx[i] - 1) > 0);
    M->val[didx] = M->idx[i + 1] - M->idx[i] - 1;
  }

  return L;
}

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
struct parilu_t *parilu_setup(uint n, const slong *const vertex, const uint nnz,
                              const uint *const row, const uint *const col,
                              const double *const val,
                              const struct parilu_opts_t *const options,
                              MPI_Comm comm, buffer *bfr) {
  struct parilu_t *ilu = tcalloc(struct parilu_t, 1);
  ilu->pivot = options->pivot;
  ilu->verbose = options->verbose;
  ilu->null_space = options->null_space;
  ilu->tol = options->tol;

  // Things to be set during factorization.
  ilu->nnz_per_row = options->nnz_per_row;
  ilu->nlvls = 0, ilu->lvl_off = NULL;
  ilu->perm = NULL;

  // Create a gslib comm out of MPI_Comm
  struct comm c;
  comm_init(&c, comm);

  // Setup CSR mat for ILU system.
  struct csr_mat *M = parilu_setup_mat(n, vertex, nnz, row, col, val, &c, bfr);
  struct csr_mat *L = parilu_setup_laplacian_mat(M);

  parilu_mat_free(M);
  parilu_mat_free(L);
  comm_free(&c);

  return ilu;
}

#undef MAX_ARRAY_SIZE
