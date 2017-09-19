/* 
 * File:   Contig_handler.h
 * Author: bx
 *
 * Created on August 26, 2017, 5:13 PM
 */

#ifndef CONTIG_HANDLER_H
#define	CONTIG_HANDLER_H

#include <stdint.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include "log.h"
#include "shc_google_sparsehash.h"
#include "encoding.h"

#define NOTHING_TO_KEEP 0
#define NOTHING_TO_REMOVE 0

// not creating vector of vector to avoid memory fragmentation
class Contig_handler {
    friend class Kmer_handler;
public:
    typedef std::vector<uint8_t>::size_type size_type;
    
    Contig_handler(char * logname);
    Contig_handler(char * logname, uint8_t rmer_length);
    
    void dump_all_contig(std::string & filename);
    
    
    // following codes change contig list 
    void push_back(uint8_t base);
    void flash();
    void delete_contig(contig_num_t i);
    void delete_contig_list(contig_num_t * remove_list, contig_num_t num);
    contig_num_t get_next_to_keep(contig_num_t * remove_list, contig_num_t remove_num, 
                                                contig_num_t remove_list_index);
    contig_num_t get_next_to_delete(contig_num_t * remove_list, contig_num_t remove_num,
                             contig_num_t keep_index, contig_num_t curr_remove_index);
    
    void declare_new_contig(kmer_count_t mean_c, size_type contig_len);
    void reject_new_contig();
    
    void print_delimitor();
    void print_delimitor(std::vector<size_type> & d);
    void print_contig_length();
    void print_contig(size_type cs, size_type len);    
    void update_delimitor_and_mean_count(contig_num_t * remove_list, contig_num_t remove_num);    
    
    void filter_contig(kmer_count_t min_count, size_type min_len, double R_threshold);
    void remove_redundance(double threshold);
    
    void get_logname(char * logname);
                
    // important contig info
    contig_num_t num_contig; 
    std::vector<uint8_t> contig_list;  
    std::vector<kmer_count_t> mean_count;        
    std::vector<size_type> delimitor; // the range is [ )   
    std::vector<size_type> contig_len_list;         
    
    //for compressing bases to byte
    uint8_t staged_num;
    uint8_t stage_count;
    
    
    //for debugging
    std::stringstream vec_str;
    char * shc_logname;            
};


#endif	/* CONTIG_HANDLER_H */

