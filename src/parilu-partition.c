#include "parilu-impl.h"
#include <math.h>

static inline void orthogonalize(scalar *const v, const uint vn,
                                 const struct comm *const c) {
  // v^T \cdot \underline{1}
  scalar sum[2] = {0, vn};
  for (uint i = 0; i < vn; i++)
    sum[0] += v[i];

  scalar wrk[2];
  comm_allreduce(c, gs_double, gs_add, sum, 2, wrk);

  // v -= (v^T \cdot \underline{1}) / ||1||^2
  sum[0] /= sum[1];
  for (uint i = 0; i < vn; i++)
    v[i] -= sum[0];
}

static inline scalar dot(const scalar *const v, const scalar *const u,
                         const uint n, const struct comm *const c) {
  // v^T \cdot u
  scalar sum = 0;
  for (uint i = 0; i < n; i++)
    sum += v[i] * u[i];

  scalar wrk;
  comm_allreduce(c, gs_double, gs_add, &sum, 1, &wrk);

  return sum;
}

static void parilu_lanczos(scalar *const fiedler, const struct parilu_mat_t *M,
                           const scalar *const init, const struct comm *const c,
                           buffer *bfr, const int verbose) {}

static void parilu_fiedler(const struct parilu_mat_t *M,
                           const struct comm *const c, buffer *bfr,
                           const int verbose) {
  parilu_debug(c, verbose, "parilu_partition: Compute Fiedler vector");

  const uint nr = M->rn;

  // Find the number of global rows and the start row id for this processor.
  slong nrg, startg;
  {
    slong out[2][1], wrk[2][1], in = nr;
    comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
    startg = out[0][0], nrg = out[1][0];
    if (nrg == 0) {
      parilu_debug(c, verbose, "parilu_partition: Number of global rows = 0");
      return;
    }
  }

  // Initialize and ortho-normalize the initial guess for Lanczos.
  scalar *init = tcalloc(scalar, nr);
  {
    for (uint i = 0; i < nr; i++)
      init[i] = startg + 1.0;

    orthogonalize(init, nr, c);

    scalar normi = 1.0 / sqrt(dot(init, init, nr, c));
    for (uint i = 0; i < nr; i++)
      init[i] *= normi;
  }
}

struct parilu_mat_t *parilu_partition(const struct parilu_mat_t *const M,
                                      const struct comm *const c, buffer *bfr,
                                      const int verbose) {
  return NULL;
}
