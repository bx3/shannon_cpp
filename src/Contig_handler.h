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

//#define LOG_CONTIG

// not creating vector of vector to avoid memory fragmentation
class Contig_handler {
    friend class Kmer_handler;
public:
    typedef std::vector<uint8_t>::size_type size_type;

    Contig_handler()
    {
        num_contig = 0;
        stage_count = 0;
        staged_num = 0x00;
        delimitor.push_back(0);
    }

    Contig_handler(bool is_com_press_)
    {
        is_use_compress = is_com_press_;
        num_contig = 0;
        stage_count = 0;
        staged_num = 0x00;
        delimitor.push_back(0);
    }

    void dump_all_contig(std::string & filename);
    void load_contig_file(std::string & filename);

    void declare_new_contig(kmer_count_t mean_c, size_type contig_len);
    void reject_new_contig();

    void print_delimitor();
    void print_contig_length();
    void print_contig(contig_num_t i);
    void log_contig(contig_num_t i, size_t len);
    void update_delimitor_and_mean_count(contig_num_t * remove_list, contig_num_t remove_num);

    void filter_contig(kmer_count_t min_count, size_type min_len, double R_threshold);
    void remove_redundance(double threshold);
    inline void get_contig(contig_num_t i, char * & base_start, uint64_t &len)
    {
        //if(is_use_compress)
        //{
        //    uint8_t * contig_start = &(contig_list.at(delimitor.at(i)));
        //    curr_contig.resize(contig_len_list.at(i));
        //    decode_byte_list(contig_start, &curr_contig.at(0), curr_contig.size());
        //    Contig_base_info contig_info(&curr_contig.at(0), curr_contig.size());
        //    return contig_info;
        //}
        //else
        //{
            base_start = (char*)(&contig_list[delimitor[i]]);
            len = contig_len_list[i];
    }

    // important contig info
    contig_num_t num_contig;
    std::vector<uint8_t> contig_list;
    std::vector<kmer_count_t> mean_count;
    std::vector<size_type> delimitor; // the range is [ )
    std::vector<size_type> contig_len_list;

    std::vector<char> curr_contig;

    //for compressing bases to byte
    uint8_t staged_num;
    uint8_t stage_count;

    bool is_use_compress;

    //for debugging
    std::stringstream vec_str;
};


#endif	/* CONTIG_HANDLER_H */
