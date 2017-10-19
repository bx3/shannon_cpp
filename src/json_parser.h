/* 
 * File:   json_parser.h
 * Author: bx
 *
 * Created on September 28, 2017, 11:02 PM
 */

#ifndef JSON_PARSER_H
#define	JSON_PARSER_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <string.h>
#include "shc_type.h"
#include "local_file_structure.h"
#include <metis.h>
#include "log.h"

struct Duplicate_setting {
    Duplicate_setting() {}
    Duplicate_setting(uint8_t rmer_length_, double threshold_, 
           kmer_count_t min_count_, size_t min_len_, bool is_use_set_):
           rmer_length(rmer_length_), threshold(threshold_), 
           min_count(min_count_), is_use_set(is_use_set_) {}
                
    void set_para(uint8_t rmer_length_, double threshold_, 
           kmer_count_t min_count_, size_t min_len_, bool is_use_set_)
    {
        rmer_length = rmer_length_;
        threshold = threshold_;
        min_count = min_count_;
        min_len = min_len_;
        is_use_set = is_use_set_;
    }
    
    uint8_t rmer_length;
    double threshold;
    kmer_count_t min_count;
    size_t min_len; 
    bool is_use_set;
};

struct Metis_setup {
    Metis_setup() {}
    Metis_setup(bool imp, idx_t partition_size_, real_t overload_, idx_t penalty_): 
            is_multiple_partition(imp), partition_size(partition_size_),  
            overload(overload_), penalty(penalty_)
    {
        inMem = 0;
        ncon = 1;
        memset(options, 0, METIS_NOPTIONS);         
    }
    //copy constructor
    Metis_setup(const Metis_setup & ms)
    {
        is_multiple_partition = ms.is_multiple_partition;
        partition_size = ms.partition_size;
        overload = ms.overload;
        penalty = ms.penalty;
        inMem = ms.inMem;
        ncon = ms.ncon;                
    }
    
    Metis_setup& operator = (const Metis_setup & ms)
    {
        is_multiple_partition = ms.is_multiple_partition;
        partition_size = ms.partition_size;
        overload = ms.overload;
        penalty = ms.penalty;
        inMem = ms.inMem;
        ncon = ms.ncon;                
    }
    
    void set_para(bool imp, idx_t partition_size_, real_t overload_, 
                                                                idx_t penalty_)            
    {
        is_multiple_partition = imp; 
        partition_size = partition_size_;
        overload = overload_; 
        penalty = penalty_;
        inMem = 0;
        ncon = 1;
        memset(options, 0, METIS_NOPTIONS);
    }
    
    idx_t partition_size;
    real_t overload;
    idx_t penalty;
    idx_t inMem;
    idx_t ncon;
    idx_t options[METIS_NOPTIONS];
    idx_t objval;
    idx_t num_vertices;
    idx_t num_partition;
    
    bool is_multiple_partition;
};

struct Output_setup {
    Output_setup(){}
    
    void set_para(bool kmer_with_info_, bool contig_, bool component_array_,
              bool read_to_component_, bool kmer_to_component_)
    {
        kmer_with_info = kmer_with_info_;
        contig = contig_;
        component_array = component_array_;
        read_to_component = read_to_component_;
        kmer_to_component = kmer_to_component_;
    }
    
    bool kmer_with_info;
    bool contig;
    bool component_array;
    bool read_to_component;
    bool kmer_to_component;
};

struct Shannon_C_setting {  
    Shannon_C_setting () {}
    
    uint8_t kmer_length;
    bool has_pair;
    bool is_compress;
    
    Local_files local_files;
    
    Duplicate_setting dup_setting;
    
    Metis_setup metis_setup;
    
    Output_setup output_setup;
};

void parser_setting_file(std::string & file_path, Shannon_C_setting & setting);
void print_setting(Shannon_C_setting & setting);
#endif	/* JSON_PARSER_H */

