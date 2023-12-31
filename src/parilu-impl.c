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

static int32_t verbose_ = 0;

/**
 * @ingroup parilu_internal_api_functions
 *
 * @brief Setup verbose level for parilu_log() function.
 *
 * @param verbose Verbosity level.
 */
void parilu_set_log_level(const int32_t verbose) { verbose_ = verbose; }

static const char *log_type_string_[3] = {"ERROR", "WARN", "INFO"};
/**
 * @ingroup parilu_internal_api_functions
 *
 * @brief Print a debug message if the verbosity level is greater than 0
 * based on the log type.
 *
 * @param c MPI communicator wrapper from gslib.
 * @param type Type of the debug message.
 * @param fmt Format string.
 * @param ... Format string arguments.
 */
void parilu_log(const struct comm *const c, const parilu_log_t type,
                const char *const fmt, ...) {
  comm_barrier(c);
  if (c->id > 0)
    return;

  if (verbose_ >= (int)type) {
    fprintf(stderr, "[%s]: ", log_type_string_[type - 1]);

    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    fflush(stderr);
    va_end(args);
  }

  if (type == PARILU_ERROR)
    MPI_Abort(c->c, EXIT_FAILURE);
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

/**
 * @ingroup parilu_user_api_functions
 *
 * @brief Finalize parILU. Frees the memory allocated for the handle returned
 * by parilu_setup() and sets the pointer to NULL.
 *
 * @param parilu Pointer to ::parilu_handle to be released.
 */
void parilu_finalize(parilu_handle **parilu) { parilu_free(parilu); }
