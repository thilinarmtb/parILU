#include "parilu-impl.h"
#include <math.h>

#define MAX_ITER 50
#define MAX_PASS 50

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

static inline void normalize(scalar *const v, const uint n,
                             const struct comm *const c) {
  scalar normi = 1.0 / sqrt(dot(v, v, n, c));
  for (uint i = 0; i < n; i++)
    v[i] *= normi;
}

static void tqli(scalar *const evec, scalar *const eval, const uint n,
                 const scalar *const alpha, const scalar *const beta) {
  return;
}

static uint lanczos_aux(scalar *const alpha, scalar *const beta,
                        scalar *const rr, const scalar *const init,
                        struct parilu_mat_op_t *op, const struct comm *const c,
                        const uint miter, buffer *bfr, const int verbose) {
  uint iter = 0;
  for (iter = 0; iter < miter; iter++) {
  }

  return iter;
}

static void parilu_lanczos(scalar *const fiedler, const struct parilu_mat_t *M,
                           const uint miter, const uint mpass,
                           const struct comm *const c, buffer *bfr,
                           const int verbose) {
  const uint rn = M->rn;
  struct parilu_mat_op_t *op = parilu_mat_op_setup(M, c, bfr);

  scalar *alpha = NULL, *beta = NULL, *rr = NULL;
  scalar *evec = NULL, *eval = NULL;
  {
    alpha = tcalloc(scalar, miter);
    beta = tcalloc(scalar, miter);
    rr = tcalloc(scalar, (miter + 1) * rn);
    evec = tcalloc(scalar, miter * miter);
    eval = tcalloc(scalar, miter);
  }

  for (uint pass = 0; pass < mpass; pass++) {
    uint iter =
        lanczos_aux(alpha, beta, rr, fiedler, op, c, iter, bfr, verbose);

    // Find eigenvalues and eigenvectors of the tridiagonal matrix.
    tqli(evec, eval, iter, alpha, beta);

    // Find min eigenvalue and associated eigenvector.
    scalar eval_min = fabs(eval[0]);
    uint eval_min_idx = 0;
    {
      for (uint i = 1; i < rn; i++) {
        if (fabs(eval[i]) < eval_min) {
          eval_min = eval[i];
          eval_min_idx = i;
        }
      }
    }

    // Project eigenvector onto the original space.
    {
      for (uint i = 0; i < rn; i++) {
        fiedler[i] = 0;
        for (uint j = 0; j < iter; j++)
          fiedler[i] += evec[j + eval_min_idx * iter] * rr[j * rn + i];
      }
      orthogonalize(fiedler, rn, c);
    }

    if (iter < miter)
      break;
  }

  // Free memory.
  {
    parilu_free(&alpha);
    parilu_free(&beta);
    parilu_free(&rr);
    parilu_free(&evec);
    parilu_free(&eval);
    parilu_mat_op_free(&op);
  }
}

static void parilu_fiedler(scalar *const fiedler, const struct parilu_mat_t *M,
                           const struct comm *const c, buffer *bfr,
                           const int verbose) {
  parilu_debug(c, verbose, PARILU_INFO,
               "parilu_partition: Compute Fiedler vector.");

  const uint nr = M->rn;

  // Find the number of global rows and the start row id for this processor.
  slong nrg, startg;
  {
    slong out[2][1], wrk[2][1], in = nr;
    comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
    startg = out[0][0], nrg = out[1][0];
    parilu_debug(c, verbose, PARILU_INFO,
                 "parilu_partition: Number of global rows = %ld.", nrg);
  }

  // Return if the global matrix is empty.
  if (nrg == 0) {
    parilu_debug(c, verbose, PARILU_WARN,
                 "parilu_partition: Number of global rows = 0.");
    return;
  }

  // Setup miter, mpass. These should be set by user.
  uint miter = MAX_ITER, mpass = MAX_PASS;
  {
    if (nrg < miter)
      miter = nrg;
    parilu_debug(c, verbose, PARILU_INFO,
                 "parilu_partition: Number of iterations = %d.", miter);
    parilu_debug(c, verbose, PARILU_INFO,
                 "parilu_partition: Number of passes = %d.", mpass);
  }

  // Initialize and ortho-normalize the initial guess for Lanczos.
  {
    for (uint i = 0; i < nr; i++)
      fiedler[i] = startg + 1.0;
    orthogonalize(fiedler, nr, c);
    normalize(fiedler, nr, c);
  }

  parilu_lanczos(fiedler, M, miter, mpass, c, bfr, verbose);

  normalize(fiedler, nr, c);
}

struct parilu_mat_t *parilu_partition(const struct parilu_mat_t *const M,
                                      const struct comm *const c, buffer *bfr,
                                      const int verbose) {
  return NULL;
}

#undef MAX_ITER
#undef MAX_PASS
