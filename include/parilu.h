#if !defined(__LIBPARILU_H__)
#define __LIBPARILU_H__

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

struct parilu_t;

/**
 * @defgroup parilu_user_api_functions User API functions
 *
 * @brief User API functions defined in `parilu.h`.
 */

PARILU_EXTERN struct parilu_t *parilu_init(int *argc, char **argv[]);

PARILU_EXTERN void parilu_finalize(struct parilu_t **parilu);

#endif // __LIBPARILU_H__
