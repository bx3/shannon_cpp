/*
 * File:   local_file_structure.h
 * Author: bx
 *
 * Created on September 18, 2017, 1:10 PM
 */

#ifndef LOCAL_FILE_STRUCTURE_H
#define	LOCAL_FILE_STRUCTURE_H
#include <iostream>
#include <syscall.h>
#include <unistd.h>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include "log.h"
#include <stdio.h>

struct Local_files;

void add_or_overwrite_directory(std::string & dir);
void add_directory_if_not_exist(std::string & dir);
void add_directory(boost::filesystem::path & dir_path);
void empty_directory(boost::filesystem::path & dir_path);
void remove_directory(boost::filesystem::path & dir_path);
void remove_file(const boost::filesystem::path & file_path);
void print_local_file_system(Local_files *lf);
void log_local_file_system(Local_files *lf);
bool exist_path(std::string a_path);
bool is_file_empty(std::string a_path);
void overwirte_a_file(std::string file_path);

size_t get_filesize(const std::string & filename);
size_t estimate_num_read(const std::string & filename, size_t read_length);
size_t estimate_num_kmer(const std::string & filename, uint8_t kmer_length);
int count_num_files(std::string path);
void copy_file(std::string file_path, std::string copy_file_path);


void convert_relative_path_to_abs(std::string rel_path, std::string & abs_path);
void convert_relative_path_to_abs(std::string & path);
bool is_abs_path(std::string & a_path);

uint64_t get_num_seq(std::string in_file);
uint64_t get_num_kmer(std::string in_file);


#define TEST_NUM_READ 1000
#define TEST_NUM_LINE 5000
#define SIZE_MULTIPLIER 1.15

struct Local_files {    //"/test_data"
    typedef std::string s;

    Local_files() {}
    Local_files(
        bool has_single_, bool has_pair_,
        s &base_path_, s &output_path_,
        s &input_kmer_path_, s &input_read_path_,
        s &input_read_path_1_, s &input_read_path_2_, s & reference_seq_path_,
        s &input_jf_path_) :
        has_single(has_single_), has_pair(has_pair_),
        base_path(base_path_), output_path(output_path_),
        input_kmer_path(input_kmer_path_), input_read_path(input_read_path_),
        input_read_path_1(input_read_path_1_), input_read_path_2(input_read_path_2_),
        reference_seq_path(reference_seq_path_),
        input_jf_path(input_jf_path_)
    {
        if(!output_path.empty())
        {
            add_output_path(output_path);
        }
    }

    void reset_paths()
    {
        if(!input_kmer_path.empty())
            has_single = true;
        else
            has_single = false;

        if(!input_read_path_1.empty() && !input_read_path_2.empty())
            has_pair = true;
        else
            has_pair = false;

        log_filename_path = output_path + "/log_shannonC";

        memcpy(shc_logname, log_filename_path.c_str(),
                                        log_filename_path.size());

        output_contig_path = output_path + "/contig";
        output_comp_path = output_path + "/comp_array";
        output_kmer_path = output_path + "/kmer";
        deleted_contig_path = output_path + "/deleted_kmer";


        output_components_read_dir = output_path + "/components_reads";
        output_components_kmer_dir = output_path + "/components_kmer";

        // Multibridge
        graph_coll_dir_name = "/comp_graph";
        output_seq_graph_path = output_path + graph_coll_dir_name;
        output_seq_graph_result_path = output_path + "/comp_graph_output";

        node_prefix = "/Graph_nodes_";
        edge_prefix = "/Graph_edges_";
        path_prefix = "/Graph_paths_";
        //SF output
        reconstructed_seq_path = output_path + "/reconstructed_seq.fasta";
        reconstructed_single_path = reconstructed_seq_path + "_single";
        reconstructed_sf_path = reconstructed_seq_path + "_sf";

        // others
        duplicate_removed_read_dir = output_path + "/duplicate_remove_reads";
        comp_read_prefix = "/comp";
        comp_kmer_prefix = "/kmer";

        eval_path = output_path + "/eval.log";
        algo_input = output_path + "/algo_input";

        unfilter_file = output_kmer_path + std::string("_unfiltered");
        sorted_unfilter_file = output_kmer_path + std::string("_unfiltered_sorted");

        single_node_dir = output_seq_graph_path + "/Single_nodes";


        // prepare dir and files
        boost::filesystem::path output_path_boost(output_path);
        if( !boost::filesystem::exists(output_path_boost) )
            add_directory(output_path_boost);
        boost::filesystem::path output_seq_graph_path_boost(output_seq_graph_path);
        if( !boost::filesystem::exists(output_seq_graph_path_boost) )
            add_directory(output_seq_graph_path_boost);
        boost::filesystem::path output_seq_graph_result_path_boost(output_seq_graph_result_path);
        //std::cout << "output_seq_graph_result_path " << output_seq_graph_result_path << std::endl;
        if( !boost::filesystem::exists(output_seq_graph_result_path_boost) )
            add_directory(output_seq_graph_result_path_boost);

        add_directory_if_not_exist(algo_input);
    }

    void add_output_path(std::string output_path_)
    {
        output_path = output_path_;
        log_filename_path = output_path + "/log_shannonC";
        output_contig_path = output_path + "/contig";
        output_comp_path = output_path + "/comp_array";
        output_kmer_path = output_path + "/kmer";
        deleted_contig_path = output_path + "/deleted_kmer";

        memcpy(shc_logname, log_filename_path.c_str(),
                                        log_filename_path.size());

        output_components_read_dir = output_path + "/components_reads";
        output_components_kmer_dir = output_path + "/components_kmer";
        // Multibridge
        graph_coll_dir_name = "/comp_graph";
        output_seq_graph_path = output_path + graph_coll_dir_name;
        output_seq_graph_result_path = output_path + "/comp_graph_output";
        node_prefix = "/Graph_nodes_";
        edge_prefix = "/Graph_edges_";
        path_prefix = "/Graph_paths_";
        //SF output
        reconstructed_seq_path = output_path + "/reconstructed_seq.fasta";
        reconstructed_single_path = reconstructed_seq_path + "_single";
        reconstructed_sf_path = reconstructed_seq_path + "_sf";

        // others
        duplicate_removed_read_dir = output_path + "/duplicate_remove_reads";
        comp_read_prefix = "/comp";
        comp_kmer_prefix = "/kmer";

        eval_path = output_path + "/eval.log";
        algo_input = output_path + "/algo_input";

        unfilter_file = output_kmer_path + std::string("_unfiltered");
        sorted_unfilter_file = output_kmer_path + std::string("_unfiltered_sorted");

        single_node_dir = output_seq_graph_path + "/Single_nodes";

        // prepare dir and files
        boost::filesystem::path output_path_boost(output_path);
        if( !boost::filesystem::exists(output_path_boost) )
            add_directory(output_path_boost);
        boost::filesystem::path output_seq_graph_path_boost(output_seq_graph_path);
        if( !boost::filesystem::exists(output_seq_graph_path_boost) )
            add_directory(output_seq_graph_path_boost);
        boost::filesystem::path output_seq_graph_result_path_boost(output_seq_graph_result_path);
        if( !boost::filesystem::exists(output_seq_graph_result_path_boost) )
            add_directory(output_seq_graph_result_path_boost);


        add_directory_if_not_exist(algo_input);
    }

    // dir path
    std::string log_filename_path;
    std::string output_comp_path;
    std::string output_path;
    std::string base_path;

    std::string reconstructed_single_path;
    std::string reconstructed_sf_path;
    //for single
    bool has_single;
    std::string input_jf_path;
    std::string input_kmer_path;
    std::string input_read_path;
    //in case paired read
    bool has_pair;
    std::string input_read_path_1;
    std::string input_read_path_2;
    // contig graph
    std::string output_contig_path;
    std::string output_kmer_path;
    std::string deleted_contig_path;
    std::string duplicate_removed_read_dir;
    // after contig graph dump
    std::string output_components_read_dir;
    std::string output_components_kmer_dir;

    std::string comp_read_prefix;
    std::string comp_kmer_prefix;
    // for dumping multibridged graph
    std::string output_seq_graph_path;
    std::string output_seq_graph_result_path;
    std::string node_prefix;
    std::string edge_prefix;
    std::string path_prefix;
    std::string single_node_dir;
    // for sparse flow
    std::string graph_coll_dir_name;
    std::string reconstructed_seq_path;

    std::string reference_seq_path;
    std::string eval_path;
    std::string algo_input;

    std::string sorted_unfilter_file;
    std::string unfilter_file;
private:
    //input name
    std::string contig_name;
    std::string comp_name;
    std::string kmer_name;
};

#endif	/* LOCAL_FILE_STRUCTURE_H */
