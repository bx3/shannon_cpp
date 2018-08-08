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
    Duplicate_setting()
    {
        rmer_length = 15;
        load_factor = 0.8;
        num_sort_thread = 1;
        threshold = 0.5;
        min_count = 3;
        min_len = 75;
        is_use_set = false;
    }
    Duplicate_setting(uint8_t rmer_length_, double threshold_,
           kmer_count_t min_count_, size_t min_len_, bool is_use_set_,
           double load_factor_, int num_sort_thread_):
           rmer_length(rmer_length_), threshold(threshold_),
           min_count(min_count_), is_use_set(is_use_set_),
           load_factor(load_factor_),
           num_sort_thread(num_sort_thread_) {}

    void set_para(uint8_t rmer_length_, double threshold_,
           kmer_count_t min_count_, size_t min_len_, bool is_use_set_,
           double load_factor_, int num_sort_thread_)
    {
        rmer_length = rmer_length_;
        threshold = threshold_;
        min_count = min_count_;
        min_len = min_len_;
        is_use_set = is_use_set_;
        load_factor = load_factor_;
        num_sort_thread = num_sort_thread_;
    }

    uint8_t rmer_length;
    double threshold;
    kmer_count_t min_count;
    size_t min_len;
    bool is_use_set;
    double load_factor;
    int num_sort_thread;
};

struct Contig_graph_setup {
    Contig_graph_setup ()
    {
        num_test = 3;
        is_assign_best = true;
        read_sampler_k = 0;
    }
    Contig_graph_setup(int num_test_, bool is_assign_best_, int read_sampler_k_) :
                       num_test(num_test_), is_assign_best(is_assign_best_),
                       read_sampler_k(read_sampler_k_) {}
    void set_para(int num_test_, bool is_assign_best_, int read_sampler_k_)
    {
        num_test=num_test_;
        is_assign_best = is_assign_best_;
        read_sampler_k = read_sampler_k_;
    }

    int num_test;
    bool is_assign_best;
    int read_sampler_k;
};

struct Metis_setup {
    Metis_setup()
    {
        is_multiple_partition = true;
        partition_size = 500;
        non_partition_size = 500;
        penalty = 5;
        overload = 2;
        inMem = 0;
        ncon = 1;
        memset(options, 0, METIS_NOPTIONS);
    }
    Metis_setup(bool imp, idx_t partition_size_, real_t overload_, idx_t penalty_,
                idx_t non_partition_size_):
            is_multiple_partition(imp), partition_size(partition_size_),
            overload(overload_), penalty(penalty_),
            non_partition_size(non_partition_size_)
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
        non_partition_size = ms.non_partition_size;
    }

    Metis_setup& operator = (const Metis_setup & ms)
    {
        is_multiple_partition = ms.is_multiple_partition;
        partition_size = ms.partition_size;
        overload = ms.overload;
        penalty = ms.penalty;
        inMem = ms.inMem;
        ncon = ms.ncon;
        non_partition_size = ms.non_partition_size;
    }

    void set_para(bool imp, idx_t partition_size_, real_t overload_,
                            idx_t penalty_, idx_t non_partition_size_)
    {
        is_multiple_partition = imp;
        partition_size = partition_size_;
        non_partition_size = non_partition_size_;
        overload = overload_;
        penalty = penalty_;
        inMem = 0;
        ncon = 1;
        memset(options, 0, METIS_NOPTIONS);
    }

    idx_t partition_size;
    idx_t non_partition_size;
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

struct Seq_graph_setup {
    Seq_graph_setup ()
    {
        max_hop_pair = 0;
        max_hop_path = 30;
        mate_pair_len = 300;
    }
    Seq_graph_setup(int max_hop_pair_,
                      int max_hop_path_,
                      int mate_pair_len_) :
                    max_hop_pair(max_hop_pair_),
                    max_hop_path(max_hop_path_),
                    mate_pair_len(mate_pair_len_) {}

    void set_para(int max_hop_pair_,
                  int max_hop_path_,
                  int mate_pair_len_)
    {
        max_hop_pair = max_hop_pair_;
        max_hop_path = max_hop_path_;
        mate_pair_len = mate_pair_len_;
    }

    int max_hop_pair;
    int max_hop_path;
    int mate_pair_len;
};

struct Sparse_flow_setup {
    Sparse_flow_setup ()
    {
        multiple_test = 8;
    }
    Sparse_flow_setup (int multiple_test_):
            multiple_test(multiple_test_) {}
    void set_para(int multiple_test_)
    {
        multiple_test = multiple_test_;
    }
    int multiple_test;
};

struct Shannon_C_setting {
    Shannon_C_setting () :
            kmer_length(0), single_read_length(0), pair_1_read_length(0),
            pair_2_read_length(0), has_single(false), has_pair(false),
            is_double_stranded(true), is_compress(false), num_parallel(1),
            take_single_node_seq(true), take_contig_seq(true){}

    bool has_single;
    uint8_t kmer_length;
    read_length_t single_read_length;

    bool has_pair;
    read_length_t pair_1_read_length;
    read_length_t pair_2_read_length;

    bool is_double_stranded;
    bool is_compress;

    //bool filter_single;
    //bool filter_paired;
    bool take_single_node_seq;
    bool take_contig_seq;

    int num_parallel;

    int output_seq_min_len;

    Local_files local_files;

    Duplicate_setting dup_setting;
    // these should have been combined
    Contig_graph_setup contig_graph_setup;
    Metis_setup metis_setup;

    Output_setup output_setup;

    Seq_graph_setup seq_graph_setup;

    Sparse_flow_setup sparse_flow_setup;
};



std::string get_setting_string(Shannon_C_setting & setting);
void parser_setting_file(std::string & file_path, Shannon_C_setting & setting);
void print_and_log_all_setting(Shannon_C_setting & setting);
void print_setting(Shannon_C_setting & setting);
void log_setting(Shannon_C_setting & setting);
void print_and_log_general_setting(Shannon_C_setting & setting);
void print_and_log_partition_setting(Shannon_C_setting & setting);
void print_and_log_multi_graph_setting(Shannon_C_setting & setting);
void print_and_log_mb_setting(Shannon_C_setting & setting);
void print_and_log_sf_setting(Shannon_C_setting & setting);
void print_and_log_general_setting(Shannon_C_setting & setting);
void print_and_log_partition_setting(Shannon_C_setting & setting);
void print_and_log_general_setting(Shannon_C_setting & setting);
void print_and_log_multi_graph_setting(Shannon_C_setting & setting);
void print_and_log_mb_setting(Shannon_C_setting & setting);
void print_and_log_sf_setting(Shannon_C_setting & setting);
void print_and_log_kmer_strand_setting(Shannon_C_setting & setting);
void print_and_log_input_path_setting(Shannon_C_setting & setting);
void print_and_log_read_length_setting(Shannon_C_setting & setting);
void print_and_log_kmer_strand_setting(Shannon_C_setting & setting);
void print_and_log_output_setting(Shannon_C_setting & setting);




#endif	/* JSON_PARSER_H */
