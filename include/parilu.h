#if !defined(__LIBPARILU_H__)
#define __LIBPARILU_H__

#include <gslib.h>

/**
 * @defgroup parilu_user_api_macros User API macros
 *
 * @brief User API macros defined in `parilu.h`.
 */

/**
 * @ingroup parilu_user_api_macros
 *
 * @def PARILU_VISIBILITY
 *
 * @brief Set the visibility of a symbol.
 * @param mode Visibility mode (hidden, default, etc.).
 */
#define PARILU_VISIBILITY(mode) __attribute__((visibility(#mode)))

/**
 * @ingroup parilu_user_api_macros
 *
 * @def PARILU_EXTERN
 *
 * @brief Declare a symbol as external.
 */
#if defined(__cplusplus)
#define PARILU_EXTERN extern "C" PARILU_VISIBILITY(default)
#else
#define PARILU_EXTERN extern PARILU_VISIBILITY(default)
#endif

/**
 * @defgroup parilu_user_api_functions User API functions
 *
 * @brief User API functions defined in `parilu.h`.
 */

/**
 * Struct to hold the options for the ILU preconditioner.
 */
struct parilu_opts_t {
  unsigned type;       /**< ILU type: ILU(0), ILUC, etc. */
  unsigned verbose;    /**< Verbosity level: 0, 1, 2, ... */
  unsigned pivot;      /**< Use pivoting or not. */
  unsigned null_space; /**< Is there a null space? */
  double
      tol; /**< 1st dropping rule: An entry a_ij is dropped abs(a_ij) < tol. */
  unsigned int nnz_per_row; /**< 2nd dropping rule: Entries are dropped so that
                               total nnz per row/col < p. */
};

PARILU_EXTERN struct parilu_opts_t *parilu_init(int *argc, char **argv[]);

struct parilu_t;

PARILU_EXTERN struct parilu_t *
parilu_setup(uint n, const slong *const vertex, const uint nnz,
             const uint *const row, const uint *const col,
             const double *const val, const struct parilu_opts_t *const options,
             MPI_Comm comm, buffer *bfr);

PARILU_EXTERN void parilu_solve(double *x, const struct parilu_t *ilu,
                                const double *b, buffer *bfr);

PARILU_EXTERN void parilu_finalize(struct parilu_t **parilu);

#endif // __LIBPARILU_H__
