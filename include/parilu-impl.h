#if !defined(__LIBPARILU_IMPL_H__)
#define __LIBPARILU_IMPL_H__

#define _POSIX_C_SOURCE 200809L

#include <parilu.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
 * parilu handle returned by parilu_init().
 */
struct parilu_t {
  unsigned verbose; /**< Verbosity level. */
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

PARILU_INTERN void parilu_debug(int verbose, const char *fmt, ...);

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

PARILU_INTERN void parilu_error(const char *fmt, ...);

#endif // __LIBPARILU_IMPL_H__
