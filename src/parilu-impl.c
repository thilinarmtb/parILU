#include "parilu-impl.h"

/**
 * @ingroup parilu_internal_api_functions
 *
 * @brief Function for freeing memory. parilu_free_() calls C standard library
 * function `free()` and sets the pointer to NULL.
 *
 * @param ptr Pointer to the memory address (pointer to a pointer) to be freed.
 */
void parilu_free_(void **ptr) { free(*ptr), *ptr = NULL; }

/**
 * @ingroup parilu_internal_api_functions
 *
 * @brief Print a debug message if the verbosity level is greater than 0.
 *
 * @param c MPI communicator wrapper from gslib.
 * @param verbose Verbosity level.
 * @param type Type of the debug message.
 * @param fmt Format string.
 * @param ... Format string arguments.
 */
void parilu_debug(const struct comm *const c, const int verbose,
                  const parilu_debug_t type, const char *const fmt, ...) {
  if (verbose > 0 && c->id == 0) {
    switch (type) {
    case PARILU_INFO:
      printf("[INFO]: ");
      break;
    case PARILU_WARN:
      printf("[WARN]: ");
      break;
    case PARILU_ERROR:
      printf("[ERROR]: ");
      break;
    default:
      break;
    }
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    printf("\n"), fflush(stdout);
    va_end(args);
  }
}

/**
 * @ingroup parilu_internal_api_functions
 *
 * @brief Assert a condition to be true and print a message if the assertion
 * fails. Use the macro parilu_assert() instead of using this function directly.
 *
 * @param cond Condition to be asserted.
 * @param msg Message to be printed if the assertion fails.
 * @param file Name of the file where the assertion is made.
 * @param line Line number where the assertion is made.
 */
void parilu_assert_(const int cond, const char *const msg,
                    const char *const file, const unsigned line) {
  if (!cond) {
    fprintf(stderr, "%s:%d Assertion failure: %s", file, line, msg);
    exit(EXIT_FAILURE);
  }
}
