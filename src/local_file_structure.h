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

struct Local_files {    //"/test_data"
    typedef std::string s;
    Local_files(s &base_path_arg, s &input_kmer_name, s &input_read_name,
               s &log_name, s &contig_name, s &comp_name, s &kmer_name) : 
               base_path_str(base_path_arg)
    {
        input_kmer_path = base_path_str + "/test_data" + input_kmer_name;
        input_read_path = base_path_str + "/test_data" + input_read_name; 
        log_filename_path = base_path_str + "/log" + log_name;
        contig_output_path = base_path_str + "/output" + contig_name;
        comp_output_path = base_path_str + "/output" + comp_name;
        kmer_output_path = base_path_str + "/output" + kmer_name;
    }
    
    std::string base_path_str;
    std::string input_kmer_path;
    std::string input_read_path;
    std::string log_filename_path;
    std::string contig_output_path;
    std::string comp_output_path;
    std::string kmer_output_path;
};


void add_directory(boost::filesystem::path & dir_path);
void remove_directory(boost::filesystem::path & dir_path);
void remove_file(const boost::filesystem::path & file_path);

#endif	/* LOCAL_FILE_STRUCTURE_H */

