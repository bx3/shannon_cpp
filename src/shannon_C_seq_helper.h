/*
 * File:   shannon_C_seq_helper.h
 * Author: bx
 *
 * Created on December 1, 2017, 2:45 PM
 */

#ifndef SHANNON_C_SEQ_HELPER_H
#define	SHANNON_C_SEQ_HELPER_H

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <algorithm>
#include <iterator>

#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include "log.h"
#include "json_parser.h"
#include "local_file_structure.h"
#include "shc_type.h"

struct JF_stats {
    JF_stats() {}
    JF_stats(uint64_t uniq_num_, uint64_t dist_num_, uint64_t total_num_,
            uint64_t max_count_) : uniq_num(uniq_num_), dist_num(dist_num_),
            total_num(total_num_), max_count(max_count_) {}
    uint64_t uniq_num;
    uint64_t dist_num;
    uint64_t total_num;
    uint64_t max_count;
};

void fasta_file_validator(std::string path);
void run_jellyfish(Shannon_C_setting & setting);
void run_jellyfish_for_file(Shannon_C_setting & setting, std::string & kmer_path);
void run_command(std::string cmd, bool print_cmd);

void parse_jf_info(Local_files & lf, JF_stats & jf_info);

void check_read_type_consistency(Shannon_C_setting & setting);
void modify_setting_with_command(Local_files & lf, boost::program_options::variables_map & vm);
int parse_main_command_line(int ac, char** av, Shannon_C_setting & setting);
uint64_t get_num_read(std::string file);

void eval_reconstructed_seq(Local_files & lf);

#ifdef USE_APPROX_COUNT
inline kmer_count_t encode_count(uint64_t count, double cr)
{
    if(count <= EXACT_THRESH || cr==1.0)
        return count;
    else
    {
        uint64_t out = EXACT_THRESH + (count-EXACT_THRESH)/cr;
        //std::cout << "out " << out << std::endl;
        assert(out<=KMER_MAX_COUNT);
        return out;
    }
}

inline uint64_t get_count(uint64_t encoded_count, double cr)
{
    if(encoded_count <= EXACT_THRESH || cr==1.0)
        return encoded_count;
    else
    {
        return EXACT_THRESH + (encoded_count-EXACT_THRESH)*cr;
    }
}
#endif

#ifdef USE_EXACT_COUNT
inline kmer_count_t encode_count(uint64_t count, double cr)
{
    return count;
}

inline uint64_t get_count(uint64_t encoded_count, double cr)
{
    return encoded_count;
}
#endif


#endif	/* SHANNON_C_SEQ_HELPER_H */
