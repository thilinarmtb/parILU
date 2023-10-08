#if !defined(__LIBPARILU_IMPL_H__)
#define __LIBPARILU_IMPL_H__

#define _POSIX_C_SOURCE 200809L

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gslib.h>

#include "parilu-defs.h"
#include "parilu-types.h"
#include "parilu.h"

/**
 * @defgroup parilu_internal_api_macros Internal API macros
 *
 * @brief Internal API macros defined in `parilu-impl.h`.
 */

/**
 * @defgroup parilu_internal_api_functions Internal API functions
 *
 * @brief Internal API functions defined in `parilu-impl.h`.
 */

/**
 * @ingroup parilu_internal_api_macros
 *
 * @def PARILU_INTERN
 *
 * @brief Declare a symbol as internal.
 */
#if defined(__cplusplus)
#define PARILU_INTERN extern "C" PARILU_VISIBILITY(hidden)
#else
#define PARILU_INTERN extern PARILU_VISIBILITY(hidden)
#endif

/**
 * parILU handle returned by parilu_setup().
 */
struct parilu_t {
  uint un;             /**< Number of local unknowns. */
  struct gs_data *u2r; /**< User vector to global vector mapping. */

  unsigned null_space; /**< Is there a null space? */
  unsigned pivot;      /**< Use pivoting or not. */

  scalar
      tol; /**< 1st dropping rule: An entry a_ij is dropped abs(a_ij) < tol. */
  uint nnz_per_row; /**< 2nd dropping rule: Entries are dropped so that total
                       nnz per row/col < p. */

  uint nlvls;    /**< Number of levels in the ILU factor. */
  uint *lvl_off; /**< Level offsets. */
  ulong *perm;   /**< Permutation vector in case of pivoting. */
  // struct par_mat_t A, L, U; /**< Matrix A. */

  struct crystal cr; /**< Crystal router for communication. */
};

/**
 * @ingroup parilu_internal_api_macros
 *
 * @brief Macro for allocating memory. This macro calls C standard library
 * function `calloc()` and casts the returned pointer as a pointer to type `T`.
 *
 * @param T Type of the memory to be allocated.
 * @param n Number of elements of type T to be allocated.
 *
 * @return Pointer to the allocated memory.
 */
#define parilu_calloc(T, n) ((T *)calloc(n, sizeof(T)))

PARILU_INTERN void parilu_free_(void **ptr);

/**
 * @ingroup parilu_internal_api_macros
 *
 * @brief Macro for freeing memory. This macro calls parilu_free_() function.
 *
 * @param p Pointer to the memory to be freed.
 */
#define parilu_free(p) parilu_free_((void **)p)

typedef enum { PARILU_INFO, PARILU_WARN, PARILU_ERROR } parilu_log_t;

PARILU_INTERN void parilu_log_init(int verbose);

PARILU_INTERN void parilu_log(const struct comm *c, parilu_log_t type,
                              const char *fmt, ...);

PARILU_INTERN void parilu_assert_(int cond, const char *fmt, const char *file,
                                  unsigned line);
/**
 * @ingroup parilu_internal_api_macros
 *
 * @brief Macro for asserting a condition. This macro calls parilu_assert_()
 * function.
 *
 * @param COND Condition to be asserted.
 * @param MSG Message to be printed if the condition is not met.
 */
#define parilu_assert(COND, MSG) parilu_assert_(COND, MSG, __FILE__, __LINE__)

/**
 * Structure for a sparse matrix.
 */
struct parilu_mat_t {
  uint rn;     /**< Local number of rows (local to processor). */
  ulong *row;  /**< Global row numbers (size: rn). */
  uint *off;   /**< Local offsets for each row (size: rn + 1). */
  uint *idx;   /**< Local column indices (size: off[rn + 1], range: [0, cn)). */
  uint cn;     /**< Local number of columns. */
  ulong *col;  /**< Global column numbers (size: cn). */
  scalar *val; /**< Local values (size: off[rn + 1]). */
};

PARILU_INTERN struct parilu_mat_t *
parilu_mat_setup(uint32_t nnz, const uint64_t *row, const uint64_t *col,
                 const double *val, const struct comm *c, int verbose);

PARILU_INTERN struct parilu_mat_t *
parilu_mat_laplacian_setup(const struct parilu_mat_t *M);

PARILU_INTERN void parilu_mat_dump(const char *file,
                                   const struct parilu_mat_t *M,
                                   const struct comm *c);

PARILU_INTERN void parilu_mat_free(struct parilu_mat_t **M);

/**
 * Structure to store the data for the matrix-vector product.
 */
struct parilu_mat_op_t {
  const struct parilu_mat_t *M; /**< Matrix. */
  struct gs_data *gsh;          /**< gs handle required for communication. */
  buffer bfr;                   /**< gslib buffer for gs workspace. */
  scalar *wrk; /**< Work array for the matrix-vector product. */
};

PARILU_INTERN struct parilu_mat_op_t *
parilu_mat_op_setup(const struct parilu_mat_t *M, const struct comm *c);

PARILU_INTERN void parilu_mat_op_apply(scalar *y, struct parilu_mat_op_t *op,
                                       const scalar *x);

PARILU_INTERN void parilu_mat_op_free(struct parilu_mat_op_t **op);

#endif // __LIBPARILU_IMPL_H__
