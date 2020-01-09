/*
 * This file is part of the CampusLink project:
 */
#include "basic_log.h"

#ifdef CL_LOGGING_ENABLED
#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include "pthread.h"

static cl_log_level filter_level = CL_LOG_LEVEL_INFO;//CL_LOG_LEVEL_VERBOSE;

/* this defines whether or not this is the first call to this function, used 
 * for resetting old log file, setting file name, etc. */
bool wlan_logfile_cleared = false;
char wlan_default_log_file[100] = "shannonC_log";

static pthread_mutex_t log_lock;


//const char * levelstring, const char * caller, 
void log_write(cl_log_level level, const char *caller, const char *format, ...)
{
    /* Only process this message if its level exceeds the current threshold */
    if (level >= filter_level)
    {
        va_list args;
        char s[1000];
        
        /* Write the log message */
        va_start(args, format);
        vsprintf(s, format, args);
        printf("[%s]%s", caller, s);
        
        //vfprintf(stderr, format, args);
        va_end(args);
    }
}

//const char * levelstring, const char * caller, 
void log_file_write(cl_log_level level, const char *caller, const char* filename,
        const char *format, ...)
{
    pthread_mutex_lock(&log_lock);
    if(!wlan_logfile_cleared){
        FILE* fp = fopen(filename, "w");
        if(fp != NULL){
            //time_t rawtime;
            //time ( &rawtime );
            //fprintf(fp, "WLAN Log File: %s\n", ctime (&rawtime));
            fclose(fp);
        }
        pthread_mutex_init(&log_lock, NULL);
        wlan_logfile_cleared = true;
    }

    FILE* fp = fopen(filename, "a+");
    if(fp!=NULL){
        /* Only process this message if its level exceeds the current threshold */
        if (level >= filter_level)
        {
            va_list args;
            char s[1000];

            //struct timespec rawtime;
            //clock_gettime(CLOCK_MONOTONIC, &rawtime);
            /* Write the log message */
            va_start(args, format);
            vsprintf(s, format, args);
            //fprintf(fp, "%.6f: [%s]%s", rawtime.tv_sec+rawtime.tv_nsec/1000000000.0, caller, s);

            if(level == CL_LOG_LEVEL_ERROR)
                printf("[%s]%s", caller, s);
            
            //vfprintf(stderr, format, args);
            va_end(args);
        }
        fclose(fp);
    }
    pthread_mutex_unlock(&log_lock);
}

#endif

