#include <gslib.h>
#include <inttypes.h>

#include <parilu.h>

static void read_matrix(uint32_t *const nnz, uint64_t **const row,
                        uint64_t **const col, double **const val,
                        const char *const file, const MPI_Comm comm,
                        const unsigned verbose) {
  struct comm c;
  comm_init(&c, comm);

  // Only rank 0 will read the file and then distribute the data to other ranks.
  // Ideally, we would like to have a parallel I/O but that is too complicated
  // for an example.
  uint64_t *row_ = NULL, *col_ = NULL;
  uint32_t count = 0;
  double *val_ = NULL;
  if (c.id == 0) {
    FILE *fp = fopen(file, "r");
    if (!fp) {
      fprintf(stderr, "Error: could not open file %s\n", file);
      MPI_Abort(comm, EXIT_FAILURE);
    }

    fscanf(fp, "%" SCNu32 "\n", nnz);

    row_ = calloc(*nnz, sizeof(**row));
    col_ = calloc(*nnz, sizeof(**col));
    val_ = calloc(*nnz, sizeof(**val));

    int64_t r, c;
    double v;
    for (uint32_t i = 0; i < *nnz; ++i) {
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
           count, *nnz);
    fflush(stdout);
  }

  if (*nnz == 0)
    return;

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
    const uint nrem = *nnz % c.np;
    const uint32_t n = (*nnz / c.np) + (nrem > 0);
    // Check invariant: n > 0 since *nnz > 0.
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

    for (uint32_t i = N; i < *nnz; ++i) {
      m.row = row_[i];
      m.col = col_[i];
      m.val = val_[i];
      m.p = (i - N) / (n - 1);
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

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 2 && rank == 0) {
    fprintf(stderr, "Usage: %s <matrix file> [<verbose level>]\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  const char *file = argv[1];
  unsigned verbose = 0;
  if (argc > 2)
    verbose = atoi(argv[2]);

  parilu_opts *opts = parilu_default_opts();
  parilu_set_verbose(opts, verbose);
  parilu_set_matrix(opts, file);

  uint32_t nnz;
  uint64_t *row = NULL, *col = NULL;
  double *val = NULL;
  read_matrix(&nnz, &row, &col, &val, file, MPI_COMM_WORLD, verbose);

  struct parilu_t *ilu = parilu_setup(nnz, row, col, val, opts, MPI_COMM_WORLD);

  free(row), free(col), free(val);
  parilu_finalize_opts(&opts);

  parilu_finalize(&ilu);

  MPI_Finalize();

  return 0;
}
