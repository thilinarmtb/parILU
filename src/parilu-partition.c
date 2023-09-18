#include "parilu-impl.h"
#include <math.h>

#define MAX_ITER 50
#define MAX_PASS 50
#define RTOL 1e-5

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
                 const scalar *const alpha, const scalar *const beta) {}

static uint lanczos_aux(scalar *const alpha, scalar *const beta,
                        scalar *const rr, const scalar *const f,
                        struct parilu_mat_op_t *op, const struct comm *const c,
                        const uint miter, const scalar rtol,
                        const int verbose) {
  const struct parilu_mat_t *M = op->M;
  const uint rn = M->rn;

  // Allocate memory for the Lanczos vectors.
  scalar *r = NULL, *p = NULL, *w = NULL;
  {
    r = tcalloc(scalar, rn);
    p = tcalloc(scalar, rn);
    w = tcalloc(scalar, rn);
  }

  // Initialize the Lanczos vectors. f should be orthogonalized wrt to
  // \underline{1} and normalized.
  {
    for (uint i = 0; i < rn; i++) {
      r[i] = f[i];
      rr[0 * rn + i] = r[i];
    }
  }

  scalar rtr = 1, rtz1 = 1, rtz2, pap1 = 0, pap2;
  uint iter = 0;
  for (iter = 1; iter <= miter; iter++) {
    rtz2 = rtz1, rtz1 = rtr;
    scalar beta_i = rtz1 / rtz2;
    if (iter == 1)
      beta_i = 0;

    // p = beta_i * p + r
    for (uint i = 0; i < rn; i++)
      p[i] = beta_i * p[i] + r[i];

    orthogonalize(p, rn, c);

    parilu_mat_op_apply(w, op, p);

    pap2 = pap1, pap1 = dot(p, w, rn, c);

    scalar alpha_i = rtz1 / pap1;
    for (uint i = 0; i < rn; i++)
      r[i] = r[i] - alpha_i * w[i];

    rtr = dot(r, r, rn, c);
    scalar rnorm = sqrt(rtr), rni = 1.0 / rnorm;
    for (uint i = 0; i < rn; i++)
      rr[iter * rn + i] = r[i] * rni;

    if (iter == 1) {
      alpha[0] = pap1;
    } else {
      alpha[iter - 1] = (pap1 + beta_i * beta_i * pap2) / rtz1;
      beta[iter - 2] = -beta_i * pap2 / sqrt(rtz2 * rtz1);
    }

    if (rnorm < rtol)
      break;
  }

  // Free memory.
  {
    parilu_free(&r);
    parilu_free(&p);
    parilu_free(&w);
  }

  return iter;
}

static void parilu_lanczos(scalar *const fiedler, const struct parilu_mat_t *M,
                           const uint miter, const uint mpass,
                           const scalar rtol, const struct comm *const c,
                           const int verbose) {
  parilu_debug(c, verbose, PARILU_INFO,
               "parilu_partition: Compute Fiedler vector.");

  const uint rn = M->rn;
  struct parilu_mat_op_t *op = parilu_mat_op_setup(M, c);

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
    parilu_debug(c, verbose, PARILU_INFO,
                 "parilu_partition: Lanczos, pass = %d.", pass);
    uint iter =
        lanczos_aux(alpha, beta, rr, fiedler, op, c, miter, rtol, verbose);

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

  // Setup miter, mpass and rtol. These should be set by user.
  uint miter = MAX_ITER, mpass = MAX_PASS;
  const scalar rtol = RTOL;
  {
    if (nrg < miter)
      miter = nrg;
    parilu_debug(c, verbose, PARILU_INFO,
                 "parilu_partition: Number of iterations = %d, Number of "
                 "passes = %d, Relative residual = %lf.",
                 miter, mpass, rtol);
  }

  // Initialize and ortho-normalize the initial guess for Lanczos.
  {
    for (uint i = 0; i < nr; i++)
      fiedler[i] = startg + 1.0;
    orthogonalize(fiedler, nr, c);
    normalize(fiedler, nr, c);
  }

  parilu_lanczos(fiedler, M, miter, mpass, rtol, c, verbose);

  normalize(fiedler, nr, c);
}

struct parilu_mat_t *parilu_partition(const struct parilu_mat_t *const M,
                                      const struct comm *const c, buffer *bfr,
                                      const int verbose) {
  return NULL;
}

#undef MAX_ITER
#undef MAX_PASS
#undef RTOL
