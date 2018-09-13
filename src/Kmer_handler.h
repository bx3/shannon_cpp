/*
 * File:   Kmer_handler.h
 * Author: bx
 *
 * Created on August 17, 2017, 5:57 PM
 */

#ifndef KMER_HANDLER_H
#define	KMER_HANDLER_H

#include <vector>
#include <map>
#include <string.h>
#include <stdlib.h>
#include <set>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <limits>
#include "encoding.h"
#include <algorithm>
#include <inttypes.h>
#include "Contig_handler.h"
#include "shc_google_sparsehash.h"
#include <bits/stringfwd.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <string>
#include <cmath>
#include <fstream>
#include "json_parser.h"
//#include "timsort.hpp"
#include <boost/lexical_cast.hpp>
#include "hopscotch_map.h"
#include "shc_setting.h"
#include "shannon_C_seq_helper.h"
//#include "parasort.h"
//#include <boost/sort/sort.hpp>


#include <sys/mman.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <vector>
#include <queue>
#include <cstdio>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <libgen.h> //for basename()


#define MAX_KMER_LENGTH 32
#define KMER_HOLDER_LEN 33
#define MAX_KMER_DICTLINE_LENGTH 50
#define DEFAULT_NUM 255
#define MAX_COUNT (std::pow(2,sizeof(kmer_count_t)*8)-1)
#define NO_MATCH 1
#define HAS_MATCH 0
#define PREFIX_OFFSET 4
//-1 is a relectant decision, otherwise my compiler gives a overflow warning,
//I guess in google sparsehash implementation, uint64_t is implicitly converted
//to a signed type
#define SPARSEHASH_DELETE_KEY (std::pow(2,(sizeof(uint64_t)*8)-1)-2)
#define SPARSEHASH_EMPTY_KEY (std::pow(2,(sizeof(uint64_t)*8)-2)-1)

//typedef std::pair<uint64_t, kmer_count_t> Kmer_Occurence_Pair;

#define KMER_LOAD_PROGRSS_STEP 5000000
#define COUT_BUFFER_FLUSH_SLEEP 100000
//#define LOG_KMER

struct Kmer_Occurence_Pair {
    Kmer_Occurence_Pair() {}
    Kmer_Occurence_Pair(uint64_t first_, kmer_count_t second_) :
                first(first_), second(second_) {}
    bool operator<(const Kmer_Occurence_Pair& kmer) const {
           return( second > kmer.second);
    }

    friend inline std::ostream& operator<<(std::ostream &os, const Kmer_Occurence_Pair &b)
    {
        os << b.first << " " << b.second;
        return os;
    }

    friend inline std::istream& operator>>(std::istream &is, Kmer_Occurence_Pair &b)
    {
        is >> b.first >> b.second;
        return is;
    }

    uint64_t first;
    kmer_count_t second;
};


struct Kmer_sorter {
    bool operator()(const Kmer_Occurence_Pair& kmer1, const Kmer_Occurence_Pair& kmer2) const {
    return( kmer1.second > kmer2.second);
    }
};

// comparison function for sorting by chromosome, then by start.
bool Kmer_sorter_func
(const Kmer_Occurence_Pair& kmer1, const Kmer_Occurence_Pair& kmer2);



class Kmer_handler {
    friend class Contig_handler;
    friend class Contig_graph_handler;
    friend class Sequence_graph_handler;
    friend class Multibridge_handler;


public:
    Kmer_handler(Local_files * lfp)
    {
        INIT_DICT(kmer_counter, SPARSEHASH_EMPTY_KEY, SPARSEHASH_DELETE_KEY);
        parse_jf_info(*lfp, jf_stats);
        num_kmer = jf_stats.dist_num;
#ifdef USE_APPROX_COUNT
        if(jf_stats.max_count <= UINT16_MAX)
        {
            compress_ratio = 1;
            std::string msg("LOADING: KMER dictionary can store count 0-65535, uses EXACT COUNT\n\t");
            msg += "based on Jellyfish stats\n";
            msg += "\t\tuniq_num  " + std::to_string(jf_stats.uniq_num) + "\n";
            msg += "\t\tdist_num  " + std::to_string(jf_stats.dist_num) + "\n";
            msg += "\t\ttotal_num  " + std::to_string(jf_stats.total_num) + "\n";
            msg += "\t\tmax_count  " + std::to_string(jf_stats.max_count)+ "\n";
            print_important_notice(msg);
        }
        else
        {
            compress_ratio = (jf_stats.max_count-EXACT_THRESH)/
                              static_cast<double>(UINT16_MAX-EXACT_THRESH);

            std::string msg("LOADING: KMER dictionary can store count 0-65535, uses APPROXIMATE COUNT\n\t");
            msg += "based on Jellyfish stats\n";
            msg += "\t\tuniq_num  " + std::to_string(jf_stats.uniq_num) + "\n";
            msg += "\t\tdist_num  " + std::to_string(jf_stats.dist_num) + "\n";
            msg += "\t\ttotal_num  " + std::to_string(jf_stats.total_num) + "\n";
            msg += "\t\tmax_count  " + std::to_string(jf_stats.max_count) + "\n";

            msg += "compression ratio " + std::to_string(compress_ratio) +
                   "But count below " +  std::to_string(EXACT_THRESH) +
                   "are stored with original count\n";
            print_important_notice(msg);
        }

        if(compress_ratio>50)
        {
            std::string msg("BE CAREFUL, TOO MUCH COMPRESSION MIGHT RESULT IN POOR PERFORMANCE \n\t");
            msg += "To change dictionary count storage range\n\t";
            msg += "Go to src/shc_type.h, comment MARCO USE_APPROX_COUNT\n\t";
            msg += "And uncomment MARCO USE_EXACT_COUNT, then compile by typing make\n";
        }
#else
        if( jf_stats.max_count > std::numeric_limits<uint32_t>::max())
        {
            shc_log_error("kmer contains count higher than 2^32-1, please use APPROX mode\n");
            exit(0);
        }

        compress_ratio = 1;
        std::string msg("LOADING: KMER dictionary can store count 0-4294967295, uses EXACT COUNT\n\t");
        msg += "based on Jellyfish stats\n";
        msg += "\t\tuniq_num  " + std::to_string(jf_stats.uniq_num) + "\n";
        msg += "\t\tdist_num  " + std::to_string(jf_stats.dist_num) + "\n";
        msg += "\t\ttotal_num  " + std::to_string(jf_stats.total_num) + "\n";
        msg += "\t\tmax_count  " + std::to_string(jf_stats.max_count) + "\n";
        print_important_notice(msg);
#endif

#if defined(USE_DENSE_KMER) || defined(USE_SPARSE_KMER)
        kmer_counter.set_deleted_key(SPARSEHASH_DELETE_KEY);
#endif
#ifdef USE_DENSE_KMER
        kmer_counter.set_empty_key(SPARSEHASH_EMPTY_KEY);
#endif
        lf = lfp;
        kmer_length = get_kmer_length_from_file(lf->output_kmer_path);

        num_kmer_deleted = 0;
        is_polyA_del = true;
        complex_thresh = 2;
        num_skip_kmer = 0;
    }

    Kmer_handler(struct Shannon_C_setting & setting_) :
                setting(setting_)
    {
        lf = &(setting_.local_files);
        kmer_length = setting_.kmer_length;
        dup_setting = setting.dup_setting;
        INIT_DICT(kmer_counter, SPARSEHASH_EMPTY_KEY, SPARSEHASH_DELETE_KEY);
        parse_jf_info(setting_.local_files, jf_stats);
        num_kmer = jf_stats.dist_num;
        if(setting.is_double_stranded)
        {
            num_kmer *= 2;
        }
        // give map a bit more space
        num_dict_kmer = static_cast<uint64_t>(
                          static_cast<double>(num_kmer)/dup_setting.load_factor);

#ifdef USE_APPROX_COUNT

        if(jf_stats.max_count <= UINT16_MAX)
        {
            compress_ratio = 1;
            std::string msg("LOADING: KMER dictionary can store count 0-65535, uses EXACT COUNT\n\t");
            msg += "based on Jellyfish stats\n";
            msg += "\t\tuniq_num  " + std::to_string(jf_stats.uniq_num) + "\n";
            msg += "\t\tdist_num  " + std::to_string(jf_stats.dist_num) + "\n";
            msg += "\t\ttotal_num  " + std::to_string(jf_stats.total_num) + "\n";
            msg += "\t\tmax_count  " + std::to_string(jf_stats.max_count)+ "\n";
            print_important_notice(msg);
        }
        else
        {
            compress_ratio = (jf_stats.max_count-EXACT_THRESH)/
                static_cast<double>(UINT16_MAX-EXACT_THRESH);

            std::string msg("LOADING: KMER dictionary can store count 0-65535, uses APPROXIMATE COUNT\n\t");
            msg += "based on Jellyfish stats\n";
            msg += "\t\tuniq_num  " + std::to_string(jf_stats.uniq_num) + "\n";
            msg += "\t\tdist_num  " + std::to_string(jf_stats.dist_num) + "\n";
            msg += "\t\ttotal_num  " + std::to_string(jf_stats.total_num) + "\n";
            msg += "\t\tmax_count  " + std::to_string(jf_stats.max_count) + "\n";

            msg += "compression ratio " + std::to_string(compress_ratio) +
                   "But count below " +  std::to_string(EXACT_THRESH) +
                   "are stored with original count\n";
            print_important_notice(msg);
        }

        if(compress_ratio>50)
        {
            std::string msg("BE CAREFUL, TOO MUCH COMPRESSION MIGHT RESULT IN POOR PERFORMANCE \n\t");
            msg += "To change dictionary count storage range\n\t";
            msg += "Go to src/shc_type.h, comment MARCO USE_APPROX_COUNT\n\t";
            msg += "And uncomment MARCO USE_EXACT_COUNT, then compile by typing make\n";
            print_warning_notice(msg);
        }
#else
        if( jf_stats.max_count > std::numeric_limits<uint32_t>::max())
        {
            shc_log_error("kmer contains count higher than 2^32-1, please use APPROX mode\n");
            exit(0);
        }
        compress_ratio = 1;
        std::string msg("LOADING: KMER dictionary can store count 0-4294967295, uses EXACT COUNT\n\t");
        msg += "based on Jellyfish stats\n";
        msg += "\t\tuniq_num  " + std::to_string(jf_stats.uniq_num) + "\n";
        msg += "\t\tdist_num  " + std::to_string(jf_stats.dist_num) + "\n";
        msg += "\t\ttotal_num  " + std::to_string(jf_stats.total_num) + "\n";
        msg += "\t\tmax_count  " + std::to_string(jf_stats.max_count) + "\n";
        print_important_notice(msg);
#endif

        if(static_cast<size_t>(num_dict_kmer) >= static_cast<size_t>(kmer_counter.max_size()))
        {
            std::cout << "num of kmer to reserve is " << num_dict_kmer << std::endl;
            std::cout << "kmer dict can at max hold " << kmer_counter.max_size()
                      << " kmers" << std::endl;
            exit(0);
        }

        RESERVE_NUM_KMER(kmer_counter, num_dict_kmer)

        if(dup_setting.is_use_set)
        {
            std::cout << "use set for error correction, reserve "
                      << (num_dict_kmer/8) << " kmer" << std::endl;
#if defined(USE_DENSE_KMER) || defined(USE_SPARSE_KMER)
                INIT_DICT(rmer_set, IMPOSSIBLE_RMER_NUM-200, IMPOSSIBLE_RMER_NUM-2);
#endif
            RESERVE_NUM_KMER(rmer_set, num_kmer/8);
        }
        else
        {
            std::cout << "use dict for error correction, reserve "
                      << (num_dict_kmer / 8) << " kmer" << std::endl;
            INIT_DICT(rmer_contig_map, IMPOSSIBLE_RMER_NUM-200, IMPOSSIBLE_RMER_NUM-2);
            RESERVE_NUM_KMER(rmer_contig_map, num_dict_kmer / 8)
        }

        num_kmer_deleted = 0;
        is_polyA_del = true;
        complex_thresh = 2;
        num_skip_kmer = 0;
    }

    void run_kmer_handler();
    void get_contig_handler(Contig_handler * c_handler);
    int build_dict_from_kmer_file();
    void sort_kmer_descending_count();
    void dump_and_sort_kmer_descending_count_external();
    void sort_kmer_descending_count_external();

    contig_num_t find_contig();
    contig_num_t find_contig_internal();
    contig_num_t find_contig_external();

    void dump_loaded_kmers_to_file(std::string & filename);
    void dump_kmers_with_info_to_file(std::string & filename);
    void load_kmers_with_info_from_file(std::string & filename);
    uint8_t get_kmer_length();
    uint8_t get_kmer_length_from_file(const std::string& filename);
    void deallocate_kmer_map();
    void clear_kmer_map();

    //debug
    uint8_t num_bit_info(uint8_t info);
    void traverse_kmer_count();
    void load_sorted_kmer(std::string sorted_kmer_file);
    void load_kmer_into_dict();
    void log_kmer(uint64_t kmer);

private:

    int add_kmer(uint64_t kmer, kmer_count_t count);
    void build_dict_from_kmer_file_helper(std::string  file);

    bool find_contig_helper(std::string & line_s, uint64_t count, bool is_rc);

    //sorting function
    bool find_prefix_kmer(uint64_t *kmer, Kmer_counter_map_iterator & next_iter);
    bool find_suffix_kmer(uint64_t *kmer, Kmer_counter_map_iterator & next_iter);
    //void kmer_change_suffix(uint64_t *kmer, char base);

    // contig relating function
    void deallocate_kmer_descend_list();
    void deallocate_contig_count_map();
    void deallocate_rmer_contig_map();
    void deallocate_rmer_set();

    bool decide_contig_and_build_rmer(double count, Contig_handler::size_type len);
    bool use_set_to_filter(kmer_count_t count, Contig_handler::size_type len);
    bool use_list_to_filter(kmer_count_t count, Contig_handler::size_type len);
    bool check_len_count(kmer_count_t count, Contig_handler::size_type len);

    void delete_kmer_for_contig(uint8_t * contig_start, Contig_handler::size_type contig_len);
    void restore_kmer_for_contig(uint8_t * contig_start, Contig_handler::size_type contig_len);

    // info concerning helper function
    bool is_kmer_info_has_ps(uint8_t info);
    bool is_kmer_info_has_prefix(uint8_t info);
    bool is_kmer_info_has_suffix(uint8_t info);
    void write_kmer_info(uint8_t num, bool is_prefix, Kmer_counter_map_iterator & kmer_it);
    void write_all_suffix_info(const uint64_t *kmer_num, Kmer_counter_map_iterator & curr_kmer_iter);
    void write_all_prefix_info(const uint64_t *kmer_num, Kmer_counter_map_iterator & curr_kmer_iter);
    bool is_info_ith_bit_set(uint8_t info, uint8_t i);

    bool is_low_complexity(std::string & base);

    kmer_len_t kmer_length;
    Kmer_counter_map kmer_counter;
    Contig_handler * ch;
    std::vector<Kmer_Occurence_Pair> kmer_descend_list;

    struct Duplicate_setting dup_setting;
    //for Rmer contig list implementation, all deallcoated
    Rmer_contig_map rmer_contig_map;
    //for Rmer set implementation, all deallcoated
    Rmer_set rmer_set;

    struct Shannon_C_setting setting;

    //for info
    size_t num_kmer_deleted;
    uint64_t num_kmer;
    uint64_t num_dict_kmer;

    Local_files *lf;

    // add to json later
    bool is_polyA_del;
    uint8_t complex_thresh;

    struct JF_stats jf_stats;
    double compress_ratio;

    Block_timer prog_timer;
    uint64_t init_kmer_size;
    uint64_t num_skip_kmer;
};

#endif	/* KMER_HANDLER_H */
