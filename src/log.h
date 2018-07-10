/*
 * File:   log.h
 * Author: bx
 *
 * Created on August 17, 2017, 12:13 AM
 */

#ifndef LOG_H

#define	LOG_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <stdarg.h>
#include <sstream>
#include <ostream>
#include <iterator>
#include "shc_type.h"
#include <time.h>
#include <iostream>
#include <pthread.h>
#include <unistd.h>

//#define LOG_KMER_DETAIL

extern int init_mem;

#define SHC_LOGGING_ENABLED
#define SHC_CONTIG_LOGGING_ENABLED

//#define LOG_KMER
//#define LOG_CONTIG
//#define LOG_DELETED_KMER
//#define LOG_SEQ_GRAPH
//#define LOG_LP_SUMMARY

#define SHOW_PROGRESS
//#define SHOW_CREATE_DIR
#define KMER_PROGRESS_STEP 1000000
#define CONTIG_PROGRESS_STEP 5000

#define MINUTE_PER_SEC 60
#define HOUR_PER_MINUTE 60

#define NANO_PER_SEC 1000000000
#define MICRO_PER_SEC 1000000
#define NANO_PER_MICRO 1000

#ifdef SHC_CONTIG_LOGGING_ENABLED
    void shc_contig_log_info(const char* format, ...);
#endif

struct Block_timer {
    Block_timer () : nTime(0), time_us(0) {}
    struct timespec nano_start;
    struct timespec nano_stamp;
    uint64_t nTime;
    int mem_start;
    uint64_t time_us;
};

extern char shc_logname[200];

template<class T> void print_bit(T a);
void print_byte_list(uint8_t* byte_list, uint16_t out_length);
void print_four_base(char* four_base, int length);
void get_vector_string(std::vector<contig_num_t> & contig_list, std::stringstream & vec_str);
void info_log_num_without_new_line(char *filename, size_t num);
void info_log_str_without_new_line(char *filename, const char * c);
void info_print_bit(char *filename, uint8_t bit_to_print);
void info_log_info(const char * filename, const char* format, ...);

void start_timer(struct Block_timer * bt);
void stop_timer(struct Block_timer * bt);
void stop_timer_np(struct Block_timer * bt); //np for no print
void print_timer(struct Block_timer * bt);

void print_important_notice(std::string msg);
void print_warning_notice(std::string msg);


inline void accurate_start_timer(struct Block_timer * bt)
{
    //
    clock_gettime(CLOCK_MONOTONIC, &bt->nano_start);
}

inline void accurate_stop_timer(struct Block_timer * bt)
{
    clock_gettime(CLOCK_MONOTONIC,&bt->nano_stamp);
    bt->nTime = (bt->nano_stamp.tv_sec - bt->nano_start.tv_sec)*NANO_PER_SEC +
             bt->nano_stamp.tv_nsec - bt->nano_start.tv_nsec;
}

inline void accurate_accumulate_timer(struct Block_timer * bt)
{
    clock_gettime(CLOCK_MONOTONIC,&bt->nano_stamp);
    bt->nTime = (bt->nano_stamp.tv_sec - bt->nano_start.tv_sec)*NANO_PER_SEC +
             bt->nano_stamp.tv_nsec - bt->nano_start.tv_nsec;
    bt->time_us += bt->nTime/NANO_PER_MICRO;
}

void log_stop_timer(struct Block_timer * bt);

int parseLine(char* line);
int get_proc_mem_value();

template<class T>
void print_bit(T a)
{
    T mask = 0x1;
    T selectBit = 0x0;
    T result = 0x0;
      //backward order
    for (int i = 0; i< sizeof(a)*8 ; i++)
    {
        selectBit = (((mask<<i) & a)>>i) & mask;
        result = (result << 1) | selectBit;
    }
    //print out bit
    for (int i = 0; i< sizeof(a)*8;i++)
    {
        printf("%lu", (((mask<<i) & result)>>i) & mask);   //print out bit in backward order
    }
    printf("\n");
}

/**
 * Severity levels for logging functions
 */
typedef enum {
    SHC_LOG_LEVEL_VERBOSE,  /**< Verbose level logging */
    SHC_LOG_LEVEL_DEBUG,    /**< Debug level logging */
    SHC_LOG_LEVEL_INFO,     /**< Information level logging */
    SHC_LOG_LEVEL_WARNING,  /**< Warning level logging */
    SHC_LOG_LEVEL_ERROR,    /**< Error level logging */
    SHC_LOG_LEVEL_CRITICAL, /**< Fatal error level logging */
    SHC_LOG_LEVEL_SILENT    /**< No output */
} shc_log_level;

/**
 * @defgroup LOG_MACROS Logging macros
 * @{
 */

/** Logs a verbose message. Does not include function/line information */
#define shc_log_verbose(FILENAME,...) \
    SHC_LOG_FILE_WRITE(FILENAME,SHC_LOG_LEVEL_VERBOSE, "[VERBOSE", __VA_ARGS__)

/** Logs a debug message. Does not include function/line information */
#define shc_log_debug(FILENAME,...) \
    SHC_LOG_FILE_WRITE(FILENAME,SHC_LOG_LEVEL_DEBUG, "[DEBUG", __VA_ARGS__)

/** Logs an info message. Does not include function/line information*/
#define shc_log_info(FILENAME,...) \
    SHC_LOG_FILE_WRITE(FILENAME,SHC_LOG_LEVEL_INFO, "[INFO" , __VA_ARGS__)

/** Logs a warning message. Includes function/line information */
#define shc_log_warning(...) \
    SHC_LOG_WRITE(SHC_LOG_LEVEL_WARNING, "[WARNING", __VA_ARGS__)

/** Logs an error message. Includes function/line information */
#define shc_log_error(...) \
    SHC_LOG_WRITE(SHC_LOG_LEVEL_ERROR, "[ERROR", __VA_ARGS__)

/** Logs a critical error message. Includes function/line information */
#define shc_log_critical(...) \
    SHC_LOG_WRITE(SHC_LOG_LEVEL_CRITICAL, "[CRITICAL", __VA_ARGS__)

#ifdef LOG_INCLUDE_FILE_INFO
#   define LOG_WRITE(LEVEL, LEVEL_STRING, ...) \
    do { log_write(LEVEL, LEVEL_STRING  \
                     " @ "  __FILE__ ":" _LOG_STRINGIFY_(__LINE__) "] " \
                     __VA_ARGS__); \
    } while (0)
#else
#   define SHC_LOG_WRITE(LEVEL, LEVEL_STRING, ...) \
    do { shc_log_write(LEVEL, __func__, LEVEL_STRING "] " __VA_ARGS__); } while (0)
#   define SHC_LOG_FILE_WRITE(FILENAME,LEVEL, LEVEL_STRING, ...) \
    do { shc_log_file_write(FILENAME, LEVEL, __func__, LEVEL_STRING "] " __VA_ARGS__); } while (0)
#endif

//declaration
#ifdef SHC_LOGGING_ENABLED
void shc_log_write(shc_log_level level, const char *caller, const char *format, ...);//writes to stdout
void shc_log_file_write(const char* filename,shc_log_level level, const char *caller,
        const char *format, ...); //writes to log files

#else
#define shc_log_write(level, caller, format, ...)  do {} while(0)
#define shc_log_file_write(filename, level, caller, format, ...)   do {} while(0)
#endif


template<class T>
void get_info_log(char * filename, T a)
{
    T mask = 0x1;
    T selectBit = 0x0;
    T result = 0x0;
      //backward order
    for (int i = 0; i< sizeof(a)*8 ; i++)
    {
        selectBit = (((mask<<i) & a)>>i) & mask;
        result = (result << 1) | selectBit;
    }
    //print out bit
    for (int i = 0; i< sizeof(a)*8;i++)
    {
        shc_log_info(filename, "%d", (((mask<<i) & result)>>i) & mask);   //print out bit in backward order
    }
}

#endif
