#include "parilu-impl.h"
#include <math.h>

#define MAX_LANCZOS_ITER 50
#define MAX_LANCZOS_PASS 50
#define MAX_TQLI_ITER 30
#define LANCZOS_RELATIVE_TOLERANCE 1e-5
#define TOLERANCE 1e-12

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

static int tqli(scalar *const evec, scalar *const eval, const uint n,
                const scalar *const alpha, const scalar *const beta) {
  if (n == 0)
    return 0;

  // Allocate and initialize workspace. Initialize evec to identity matrix.
  scalar *d = NULL, *e = NULL;
  {
    d = parilu_calloc(scalar, n);
    for (uint i = 0; i < n; i++)
      d[i] = alpha[i];

    e = parilu_calloc(scalar, n);
    for (uint i = 0; i < n - 1; i++)
      e[i] = beta[i];

    for (uint i = 0; i < n; i++) {
      for (uint j = 0; j < n; j++)
        evec[i * n + j] = 0;
      evec[i * n + i] = 1;
    }
  }

  int err = 0;
  uint m;
  for (uint l = 0; l < n; l++) {
    uint iter = 0;
    do {
      for (m = l; m < n - 1; m++) {
        scalar dd = fabs(d[m]) + fabs(d[m + 1]);
        if (fabs(e[m]) / dd < TOLERANCE)
          break;
      }

      if (m != l) {
        for (uint i = 0; i < n; i++)
          eval[i] = d[i];
        if (iter++ == MAX_TQLI_ITER) {
          err = 1;
          goto cleanup;
        }
      }

      scalar g = (d[l + 1] - d[l]) / (2 * e[l]);
      scalar r = sqrt(g * g + 1);
      g = d[m] - d[l] + e[l] / (g + copysign(r, g));

      scalar s = 1, c = 1, p = 0;
      sint i;
      for (i = m - 1; i >= (sint)l; i--) {
        scalar f = s * e[i], b = c * e[i];
        if (fabs(f) >= fabs(g)) {
          c = g / f;
          r = sqrt(c * c + 1);
          e[i + 1] = f * r;
          s = 1 / r;
          c *= s;
        } else {
          s = f / g;
          r = sqrt(s * s + 1);
          e[i + 1] = g * r;
          c = 1 / r;
          s *= c;
        }

        g = d[i + 1] - p;
        r = (d[i] - g) * s + 2 * c * b;
        p = s * r;
        d[i + 1] = g + p;
        g = c * r - b;

        for (uint k = 0; k < n; k++) {
          f = evec[k * n + i + 1];
          evec[k * n + i + 1] = s * evec[k * n + i] + c * f;
          evec[k * n + i] = c * evec[k * n + i] - s * f;
        }
      }

      if (r < TOLERANCE && i >= (sint)l)
        continue;

      d[l] -= p;
      e[l] = g;
      e[m] = 0;
    } while (m != l);
  }

  // Transpose the vectors in evec to match C ordering.
  for (uint i = 0; i < n; i++) {
    for (uint j = 0; j < i; j++) {
      scalar tmp = evec[i * n + j];
      evec[i * n + j] = evec[j * n + i];
      evec[j * n + i] = tmp;
    }
  }

  // Normalize eigenvectors and copy eigenvalues.
  for (uint k = 0; k < n; k++) {
    e[k] = 0;
    for (uint i = 0; i < n; i++)
      e[k] += evec[k * n + i] * evec[k * n + i];

    if (e[k] <= TOLERANCE) {
      err = 1;
      goto cleanup;
    }

    e[k] = sqrt(fabs(e[k]));
    scalar scale = 1.0 / e[k];
    for (uint i = 0; i < n; i++)
      evec[i * n + k] *= scale;
  }

  for (uint i = 0; i < n; i++)
    eval[i] = d[i];

cleanup:
  parilu_free(&e);
  parilu_free(&d);

  return err;
}

static uint lanczos_aux(scalar *const alpha, scalar *const beta,
                        scalar *const rr, const scalar *const f,
                        parilu_matrix_operator *op, const struct comm *const c,
                        const uint miter, const scalar rtol) {
  parilu_log(c, PARILU_INFO, "lanczos_aux: miter = %d, rtol = %e", miter, rtol);

  const parilu_matrix *M = op->M;
  const uint rn = M->rn;

  // Allocate memory for the Lanczos vectors.
  scalar *r = NULL, *p = NULL, *w = NULL;
  {
    r = parilu_calloc(scalar, rn);
    p = parilu_calloc(scalar, rn);
    w = parilu_calloc(scalar, rn);
  }

  // Initialize the Lanczos vectors. f should be orthogonalized wrt to
  // \underline{1} and normalized.
  for (uint i = 0; i < rn; i++) {
    r[i] = f[i];
    rr[0 * rn + i] = r[i];
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

    parilu_matrix_operator_apply(w, op, p);

    pap2 = pap1, pap1 = dot(p, w, rn, c);

    scalar alpha_i = rtz1 / pap1;
    for (uint i = 0; i < rn; i++)
      r[i] = r[i] - alpha_i * w[i];

    rtr = dot(r, r, rn, c);
    scalar rnorm = sqrt(rtr), rni = 1.0 / rnorm;
    for (uint i = 0; i < rn; i++)
      rr[iter * rn + i] = r[i] * rni;

    if (iter == 1) {
      alpha[iter - 1] = pap1 / rtz1;
    } else {
      alpha[iter - 1] = (pap1 + beta_i * beta_i * pap2) / rtz1;
      beta[iter - 2] = -beta_i * pap2 / sqrt(rtz2 * rtz1);
    }

    parilu_log(c, PARILU_INFO, "lanczos_aux: iter = %02d rnorm = %e", iter,
               rnorm);
    if (rnorm < rtol)
      break;
  }

  // Free memory.
  {
    parilu_free(&r);
    parilu_free(&p);
    parilu_free(&w);
  }

  return iter <= miter ? iter : miter;
}

static void parilu_lanczos(scalar *const fiedler, const parilu_matrix *M,
                           const uint miter, const uint mpass,
                           const scalar rtol, const struct comm *const c) {
  parilu_log(c, PARILU_INFO, "parilu_lanczos: Run Lanczos ...");

  // Find the number of global rows and the start row id for this processor.
  const uint rn = M->rn;
  slong rng, startg;
  {
    slong out[2][1], wrk[2][1], in = rn;
    comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
    startg = out[0][0], rng = out[1][0];

    // Return if the global matrix is empty.
    if (rng == 0) {
      parilu_log(c, PARILU_WARN, "parilu_lanczos: Number of global rows = 0");
      return;
    }

    // Return if the matrix is a scalar.
    parilu_log(c, PARILU_INFO, "parilu_lanczos: Number of global rows = %ld",
               rng);
    if (rng == 1) {
      if (rn == 1)
        fiedler[0] = 1.0;
      return;
    }
  }

  uint miter_ = miter;
  {
    if (rng < miter)
      miter_ = rng;
    parilu_log(c, PARILU_INFO,
               "parilu_lanczos: Number of iterations = %d, Number of "
               "passes = %d, Relative residual = %e",
               miter_, mpass, rtol);

    // Initialize the fiedler vector.
    for (uint i = 0; i < rn; i++)
      fiedler[i] = startg + i + 1;
  }

  parilu_matrix_operator *op = parilu_matrix_operator_setup(M, c->c);

  scalar *alpha = NULL, *beta = NULL, *rr = NULL;
  scalar *evec = NULL, *eval = NULL;
  {
    alpha = parilu_calloc(scalar, miter);
    beta = parilu_calloc(scalar, miter);
    rr = parilu_calloc(scalar, (miter + 1) * rn);
    evec = parilu_calloc(scalar, miter * miter);
    eval = parilu_calloc(scalar, miter);
  }

  for (uint pass = 0; pass < mpass; pass++) {
    parilu_log(c, PARILU_INFO, "parilu_partition: Lanczos, pass = %d",
               pass + 1);
    uint iter = lanczos_aux(alpha, beta, rr, fiedler, op, c, miter_, rtol);

    // Find eigenvalues and eigenvectors of the tridiagonal matrix.
    sint err = tqli(evec, eval, iter, alpha, beta), wrk;
    comm_allreduce(c, gs_int, gs_add, &err, 1, &wrk);
    if (err) {
      parilu_log(c, PARILU_ERROR, "parilu_partition: tqli failed");
      goto cleanup;
    }

    // Find min eigenvalue and associated eigenvector.
    scalar eval_min = fabs(eval[0]);
    uint eval_min_idx = 0;
    for (uint i = 1; i < iter; i++) {
      if (fabs(eval[i]) < eval_min) {
        eval_min = fabs(eval[i]);
        eval_min_idx = i;
      }
    }

    // Project eigenvector onto the original space.
    for (uint i = 0; i < rn; i++) {
      fiedler[i] = 0;
      for (uint j = 0; j < iter; j++)
        fiedler[i] += evec[j + eval_min_idx * iter] * rr[j * rn + i];
    }

    orthogonalize(fiedler, rn, c);
    normalize(fiedler, rn, c);

    if (iter < miter)
      break;
  }

cleanup:
  // Free memory.
  {
    parilu_free(&alpha);
    parilu_free(&beta);
    parilu_free(&rr);
    parilu_free(&evec);
    parilu_free(&eval);
    parilu_matrix_operator_free(&op);
  }
}

static void parilu_fiedler(const parilu_matrix *M, const struct comm *const c,
                           buffer *bfr) {
  parilu_log(c, PARILU_INFO, "parilu_partition: Compute Fiedler vector");

  int miter = MAX_LANCZOS_ITER, mpass = MAX_LANCZOS_PASS;
  scalar rtol = LANCZOS_RELATIVE_TOLERANCE;

  scalar *fiedler = parilu_calloc(scalar, M->rn);
  parilu_lanczos(fiedler, M, miter, mpass, rtol, c);
  parilu_free(&fiedler);
}

void parilu_partition(const parilu_matrix *const M, const MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  parilu_log(&c, PARILU_INFO, "parilu_partition: Partition matrix.");

  buffer bfr;
  buffer_init(&bfr, 1024);

  parilu_fiedler(M, &c, &bfr);

  buffer_free(&bfr);

  parilu_log(&c, PARILU_INFO, "parilu_partition: done.");

  comm_free(&c);
}

#undef MAX_LANCZOS_ITER
#undef MAX_LANCZOS_PASS
#undef MAX_TQLI_ITER
#undef LANCZOS_RELATIVE_TOLERANCE
#undef TOLERANCE
