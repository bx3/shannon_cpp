/* 
 * File:   basic_log.h
 * Author: bx
 *
 * Created on September 1, 2017, 1:38 PM
 */

#ifndef BASIC_LOG_H
#define	BASIC_LOG_H

#define CL_LOGGING_ENABLED
//#define CL_LOG_RX_PIPE


#define _LOG_STRINGIFY__(x)         #x
#define _LOG_STRINGIFY_(x)          _LOG_STRINGIFY__(x)

extern char shannonC_default_log_file[100];
extern bool shannonC_logfile_cleared;
/**
 * Severity levels for logging functions
 */
typedef enum {
    CL_LOG_LEVEL_VERBOSE,  /**< Verbose level logging */
    CL_LOG_LEVEL_DEBUG,    /**< Debug level logging */
    CL_LOG_LEVEL_INFO,     /**< Information level logging */
    CL_LOG_LEVEL_WARNING,  /**< Warning level logging */
    CL_LOG_LEVEL_ERROR,    /**< Error level logging */
    CL_LOG_LEVEL_CRITICAL, /**< Fatal error level logging */
    CL_LOG_LEVEL_SILENT    /**< No output */
} cl_log_level;

/**
 * @defgroup LOG_MACROS Logging macros
 * @{
 */

/** Logs a verbose message. Does not include function/line information */
#define log_verbose(...) \
    LOG_FILE_WRITE(CL_LOG_LEVEL_VERBOSE, "[VERBOSE", __VA_ARGS__)

/** Logs a debug message. Does not include function/line information */
#define log_debug(...) \
    LOG_FILE_WRITE(CL_LOG_LEVEL_DEBUG, "[DEBUG", __VA_ARGS__)

/** Logs an info message. Does not include function/line information*/
#define log_info(...) \
    LOG_FILE_WRITE(CL_LOG_LEVEL_INFO, "[INFO" , __VA_ARGS__)

/** Logs a warning message. Includes function/line information */
#define log_warning(...) \
    LOG_WRITE(CL_LOG_LEVEL_WARNING, "[WARNING", __VA_ARGS__)

/** Logs an error message. Includes function/line information */
#define log_error(...) \
    LOG_FILE_WRITE(CL_LOG_LEVEL_ERROR, "[ERROR", __VA_ARGS__)


/** Logs a critical error message. Includes function/line information */
#define log_critical(...) \
    LOG_WRITE(CL_LOG_LEVEL_CRITICAL, "[CRITICAL", __VA_ARGS__)

#ifdef CL_LOG_RX_PIPE
#define log_rx_pipe(...) \
        LOG_FILE_WRITE(CL_LOG_LEVEL_INFO, "[INFO", __VA_ARGS__)
#else
    #define log_rx_pipe(...) \
            do{} while (0)
#endif
    
            
/** @} */

#ifdef LOG_INCLUDE_FILE_INFO
#   define LOG_WRITE(LEVEL, LEVEL_STRING, ...) \
    do { log_write(LEVEL, LEVEL_STRING  \
                     " @ "  __FILE__ ":" _LOG_STRINGIFY_(__LINE__) "] " \
                     __VA_ARGS__); \
    } while (0)
#else
#   define LOG_WRITE(LEVEL, LEVEL_STRING, ...) \
    do { log_write(LEVEL, __func__, LEVEL_STRING "] " __VA_ARGS__); } while (0)
#   define LOG_FILE_WRITE(LEVEL, LEVEL_STRING, ...) \
    do { log_file_write(LEVEL, __func__, wlan_default_log_file, LEVEL_STRING "] " __VA_ARGS__); } while (0)
#endif

/**
 * Writes a message to the logger with the specified log level. This function
 * should only be invoked indirectly through one of the log_*
 * convenience macros above to ensure consistent formatting of log messages.
 *
 * @param   level       The severity level of the message
 * @param   format      The printf-style format string
 * @param   ...         Optional values for format specifier substitution
 */
#ifdef CL_LOGGING_ENABLED
//void log_write(cl_log_level level, const char *format, ...);
void log_write(cl_log_level level, const char *caller, const char *format, ...);//writes to stdout
void log_file_write(cl_log_level level, const char *caller, const char* filename, 
        const char *format, ...); //writes to log files
#else
#define log_write(level, format, ...)  do {} while(0)
#define log_file_write(level, caller, filename, format, ...)   do {} while(0)
#endif

#endif	/* BASIC_LOG_H */

