#include <inttypes.h>
#include <stdio.h>

#include "parilu-impl.h"

void parilu_read_matrix(uint32_t *const nnz, uint64_t **const row,
                        uint64_t **const col, double **const val,
                        const char *const file, const MPI_Comm comm,
                        const unsigned verbose) {
  struct comm c;
  comm_init(&c, comm);

  // Only rank 0 will read the file and then distribute the data to other ranks.
  // Ideally, we would like to have a parallel I/O but that is too complicated
  // for an example.
  uint64_t *row_ = NULL, *col_ = NULL;
  uint32_t count = 0, nnz_ = 0;
  double *val_ = NULL;
  if (c.id == 0) {
    FILE *fp = fopen(file, "r");
    if (!fp) {
      fprintf(stderr, "Error: could not open file %s\n", file);
      MPI_Abort(comm, EXIT_FAILURE);
    }

    fscanf(fp, "%" SCNu32 "\n", &nnz_);

    row_ = calloc(nnz_, sizeof(**row));
    col_ = calloc(nnz_, sizeof(**col));
    val_ = calloc(nnz_, sizeof(**val));

    int64_t r, c;
    double v;
    for (uint32_t i = 0; i < nnz_; ++i) {
      fscanf(fp, "%" SCNd64 " %" SCNd64 " %lf\n", &r, &c, &v);
      if (r < 0 || c < 0)
        continue;
      row_[count] = r, col_[count] = c, val_[count] = v;
      count++;
    }

    fclose(fp);
  }

  if (c.id == 0 && verbose) {
    printf("read_matrix: Read %" PRIu32 " nonzeros from %" PRIu32 " entries.\n",
           count, nnz_);
    fflush(stdout);
  }

  struct entry_t {
    uint64_t row;
    uint64_t col;
    double val;
    uint32_t p;
  };

  // Fill the matrix data at rank 0 to an array.
  struct array mat;
  array_init(struct entry_t, &mat, 1024);
  if (c.id == 0) {
    // Size of the expected data on rank 0.
    const uint nrem = count % c.np;
    const uint32_t n = (count / c.np) + (nrem > 0);
    // Check invariant: n > 0 since count > 0.
    assert(n >= 1);

    struct entry_t m;
    const uint32_t N = n * nrem;
    for (uint32_t i = 0; i < N; ++i) {
      m.row = row_[i];
      m.col = col_[i];
      m.val = val_[i];
      m.p = i / n;
      array_cat(struct entry_t, &mat, &m, 1);
    }

    for (uint32_t i = N; i < count; ++i) {
      m.row = row_[i];
      m.col = col_[i];
      m.val = val_[i];
      m.p = (i - N) / (n - nrem);
      array_cat(struct entry_t, &mat, &m, 1);
    }
  }
  free(row_), free(col_), free(val_);

  // Send the data to other ranks.
  {
    struct crystal cr;
    crystal_init(&cr, &c);
    sarray_transfer(struct entry_t, &mat, p, 0, &cr);
    crystal_free(&cr);
  }

  *nnz = mat.n;
  if (mat.n > 0) {
    row_ = *row = malloc(mat.n * sizeof(**row));
    col_ = *col = malloc(mat.n * sizeof(**col));
    val_ = *val = malloc(mat.n * sizeof(**val));
    const struct entry_t *const ptr = (const struct entry_t *const)mat.ptr;
    for (uint32_t i = 0; i < mat.n; ++i) {
      row_[i] = ptr[i].row;
      col_[i] = ptr[i].col;
      val_[i] = ptr[i].val;
    }
  }

  array_free(&mat), comm_free(&c);
}
