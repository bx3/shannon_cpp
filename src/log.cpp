#include "log.h"

// for debugging
void get_vector_string(std::vector<contig_num_t> & contig_list, std::stringstream & vec_str)
{        
    vec_str.str("");
    std::copy(contig_list.begin(), contig_list.end(), std::ostream_iterator<int>(vec_str, " "));    
}


void print_byte_list(uint8_t* byte_list, uint16_t out_length)
{
    for (int i=0; i< out_length; i++)
        printf("%u ", byte_list[i]);
    printf("\n");
}

void print_four_base(char* four_base, int length)
{
    char *base = new char[length+1];
    memcpy(base, four_base, length);
    base[length] = '\0';
    printf("%s \t", base);
}

//static bool info_logfile_cleared = false;

void info_log_info(char * filename, const char* format, ...)
{
    //char filename[50];       
    //if(!info_logfile_cleared){
    //    fclose(fopen(filename, "w"));
    //    info_logfile_cleared = true;
    //}      
    va_list args;
    va_start (args, format);   
    FILE* stream;
    stream = fopen(filename,"a+");     
    vfprintf (stream, format, args);
    va_end (args); 
    fclose(stream);  
}

// hard coded using with caution
void info_log_num_without_new_line(char *filename, size_t num)
{
    info_log_info(filename, "%3d ", num);            
}

void info_log_str_without_new_line(char *filename, const char * c)
{    
    info_log_info(filename, "%s ", c);            
}

void info_print_bit(char *filename, uint8_t bit_to_print)
{
    uint16_t mask = 0x1;
    uint16_t selectBit = 0x0; 
    uint16_t result = 0x0;
    //backward order
    for (int i = 0; i< 8 ; i++)
    {
        selectBit = (((mask<<i) & bit_to_print)>>i) & mask;
        result = (result << 1) | selectBit;  
    } 
    //print out bit
    info_log_info(filename, "[Info content][INFO] ");
    for (int i = 0; i< 8;i++)
    {
        info_log_info(filename, "%d", (((mask<<i) & result)>>i) & mask);
    } 
    info_log_info(filename, "\n");
}

#ifdef SHC_LOGGING_ENABLED
#include <time.h>

static shc_log_level filter_level = SHC_LOG_LEVEL_VERBOSE;
static bool logfile_cleared = false;

/** this function is used with MARCO to write the desired message on the screen*/
void shc_log_write(shc_log_level level, const char *caller, const char *format, ...)
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

/** 
 * this function is used with MARCO to write data in a certain file specified by the
 * filename in the argument. The file will be automatically cleared when this is a 
 * new process/test. 
 */
void shc_log_file_write(const char* filename,shc_log_level level, const char *caller, 
        const char *format, ...)
{
    if(!logfile_cleared){
        FILE* fp = fopen(filename, "w");
        if(fp != NULL){
            time_t rawtime;
            time ( &rawtime );
            fprintf(fp, "SHC Log File: [%s] \n", ctime (&rawtime));
            fclose(fp);
        }
        logfile_cleared = true;
    }

    FILE* fp = fopen(filename, "a+");
    if(fp!=NULL){
        /* Only process this message if its level exceeds the current threshold */
        if (level >= filter_level)
        {
            va_list args;
            char s[1000];

            /* Write the log message */
            va_start(args, format);
            vsprintf(s, format, args);
            fprintf(fp, "[%s]%s", caller, s);
            va_end(args);
        }
        fclose(fp);
    }
}


#ifdef DCF_TX_CONTENT_LOGGING_ENABLED
/**
 * this function is used to write detailed  information into a log
 * file, which includes delimitor
 */
#include <time.h>
static bool sch_contig_logfile_cleared = false;
void tx_content_log_info(uint8_t* addr, const char* format, ...)
{
    char filename[50];   
    sprintf(filename, "TX_CONTENT_%02x%02x%02x%02x%02x%02x.log", addr[0], addr[1], addr[2], addr[3], addr[4], addr[5]);
    if(!tx_content_logfile_cleared){
        fclose(fopen(filename, "w"));
        tx_content_logfile_cleared = true;
    }      
    va_list args;
    va_start (args, format);   
    FILE* stream;
    stream = fopen(filename,"a+");     
    vfprintf (stream, format, args);
    va_end (args); 
    fclose(stream);  
}   

/**
 * this function is used with tx_content_log_info, this function is used to 
 * print out every bit of a uint16_t type, the LSB is on the right side.
 * @param addr
 * @param bit_to_print
 */
void tx_content_print_bit(uint8_t *addr, uint16_t bit_to_print)
{
    uint16_t mask = 0x1;
    uint16_t selectBit = 0x0; 
    uint16_t result = 0x0;
    //backward order
    for (int i = 0; i< 16 ; i++)
    {
        selectBit = (((mask<<i) & bit_to_print)>>i) & mask;
        result = (result << 1) | selectBit;  
    } 
    //print out bit
    for (int i = 0; i< 16;i++)
    {
        tx_content_log_info(addr, "%d", (((mask<<i) & result)>>i) & mask);
    } 
}
#endif

#ifdef SHC_CONTIG_LOGGING_ENABLED
/**
 * this function is used to write detailed  information into a log
 * file, which includes delimitor
 */
static bool shc_contig_logfile_cleared = false;
void shc_contig_log_info(const char* format, ...)
{
    char filename[50];   
    sprintf(filename, "SHC_CONIG_INFO.log");
    if(!shc_contig_logfile_cleared){
        fclose(fopen(filename, "w"));
        shc_contig_logfile_cleared = true;
    }      
    va_list args;
    va_start (args, format);   
    FILE* stream;
    stream = fopen(filename,"a+");     
    vfprintf (stream, format, args);
    va_end (args); 
    fclose(stream);  
}   
#endif

#endif

