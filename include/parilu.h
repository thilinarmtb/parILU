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

PARILU_EXTERN void parilu_read_matrix(uint32_t *nnz, uint64_t **row,
                                      uint64_t **col, double **val,
                                      const char *file, MPI_Comm comm,
                                      unsigned verbose);
/**
 * @typedef parilu_options
 *
 * @brief Store the options for parILU. This is a typedef for struct
 * ::parilu_options_t.
 */
typedef struct parilu_options_t parilu_options;
PARILU_EXTERN parilu_options *parilu_default_options(void);

PARILU_EXTERN int parilu_set_verbose(parilu_options *options, unsigned verbose);
PARILU_EXTERN int parilu_get_verbose(const parilu_options *options,
                                     unsigned *verbose);

PARILU_EXTERN int parilu_set_type(parilu_options *options, unsigned type);
PARILU_EXTERN int parilu_get_type(const parilu_options *options,
                                  unsigned *type);

PARILU_EXTERN int parilu_set_pivot(parilu_options *options, unsigned pivot);
PARILU_EXTERN int parilu_get_pivot(const parilu_options *options,
                                   unsigned *pivot);

PARILU_EXTERN int parilu_set_null_space(parilu_options *options,
                                        unsigned null_space);
PARILU_EXTERN int parilu_get_null_space(const parilu_options *options,
                                        unsigned *null_space);

PARILU_EXTERN int parilu_set_tol(parilu_options *options, double tol);
PARILU_EXTERN int parilu_get_tol(const parilu_options *options, double *tol);

PARILU_EXTERN int parilu_set_nnz_per_row(parilu_options *options,
                                         unsigned nnz_per_row);
PARILU_EXTERN int parilu_get_nnz_per_row(const parilu_options *options,
                                         unsigned *nnz_per_row);

PARILU_EXTERN void parilu_finalize_options(parilu_options **options);

/**
 * @typedef parilu_handle
 *
 * @brief parILU handle returned by parilu_setup(). This is a typedef of struct
 * ::parilu_handle_t.
 */
typedef struct parilu_handle_t parilu_handle;
PARILU_EXTERN parilu_handle *
parilu_setup(uint32_t nnz, const uint64_t *row, const uint64_t *col,
             const double *val, const parilu_options *options, MPI_Comm comm);

PARILU_EXTERN void parilu_solve(double *x, const parilu_handle *ilu,
                                const double *b);

PARILU_EXTERN void parilu_finalize(parilu_handle **parilu);

#endif // __LIBPARILU_H__
