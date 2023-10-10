#if !defined(__LIBPARILU_H__)
#define __LIBPARILU_H__

#include <mpi.h>
#include <stdint.h>

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
 * Store the options for parILU. This is a typedef for struct ::parilu_opts_t.
 */
typedef struct parilu_opts_t parilu_opts;
PARILU_EXTERN parilu_opts *parilu_default_opts(void);

PARILU_EXTERN int parilu_set_verbose(parilu_opts *opts, unsigned verbose);
PARILU_EXTERN int parilu_get_verbose(const parilu_opts *opts,
                                     unsigned *verbose);

PARILU_EXTERN int parilu_set_type(parilu_opts *opts, unsigned type);
PARILU_EXTERN int parilu_get_type(const parilu_opts *opts, unsigned *type);

PARILU_EXTERN int parilu_set_pivot(parilu_opts *opts, unsigned pivot);
PARILU_EXTERN int parilu_get_pivot(const parilu_opts *opts, unsigned *pivot);

PARILU_EXTERN int parilu_set_null_space(parilu_opts *opts, unsigned null_space);
PARILU_EXTERN int parilu_get_null_space(const parilu_opts *opts,
                                        unsigned *null_space);

PARILU_EXTERN int parilu_set_tol(parilu_opts *opts, double tol);
PARILU_EXTERN int parilu_get_tol(const parilu_opts *opts, double *tol);

PARILU_EXTERN int parilu_set_nnz_per_row(parilu_opts *opts,
                                         unsigned nnz_per_row);
PARILU_EXTERN int parilu_get_nnz_per_row(const parilu_opts *opts,
                                         unsigned *nnz_per_row);

PARILU_EXTERN int parilu_set_matrix(parilu_opts *opts, const char *file);
PARILU_EXTERN int parilu_get_matrix(const parilu_opts *opts, char **file);

PARILU_EXTERN void parilu_finalize_opts(parilu_opts **opts);

/**
 * parILU handle returned by parilu_setup(). This is a typedef of struct
 * ::parilu_t.
 */
typedef struct parilu_t parilu;
PARILU_EXTERN parilu *parilu_setup(uint32_t nnz, const uint64_t *row,
                                   const uint64_t *col, const double *val,
                                   const parilu_opts *options, MPI_Comm comm);

PARILU_EXTERN void parilu_solve(double *x, const parilu *ilu, const double *b);

PARILU_EXTERN void parilu_finalize(parilu **parilu);

#endif // __LIBPARILU_H__
