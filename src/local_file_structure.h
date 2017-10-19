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
#include <fstream>

struct Local_files;

void add_directory(boost::filesystem::path & dir_path);
void replace_directory(boost::filesystem::path & dir_path);
void remove_directory(boost::filesystem::path & dir_path);
void remove_file(const boost::filesystem::path & file_path);
void print_local_file_system(Local_files *lf);
size_t get_filesize(const std::string & filename);
size_t estimate_num_read(const std::string & filename, size_t read_length);
size_t estimate_num_kmer(const std::string & filename, uint8_t kmer_length);

#define TEST_NUM_READ 1000
#define TEST_NUM_LINE 5000
#define SIZE_MULTIPLIER 1.15

struct Local_files {    //"/test_data"
    typedef std::string s;
    
    Local_files() {}
    
    // for single stranded read
    Local_files(s &base_path_arg, s &output_dir_name, s &input_kmer_name, s &input_read_name,
               s &log_name, s &contig_name_arg, s &comp_name_arg, s &kmer_name_arg) : 
               base_path_str(base_path_arg) , contig_name(contig_name_arg), 
               comp_name(comp_name_arg), kmer_name(kmer_name_arg)
    {
        output_path_str = base_path_str + output_dir_name;
        input_kmer_path = base_path_str + "/test_data" + input_kmer_name;
        input_read_path = base_path_str + "/test_data" + input_read_name; 
        log_filename_path = base_path_str + "/log" + log_name;
        output_contig_path = output_path_str + contig_name;
        output_comp_path = output_path_str + comp_name;
        output_kmer_path = output_path_str + kmer_name;
        deleted_contig_path = output_path_str + "/deleted_kmer";  
                    
        output_components_read_dir = output_path_str + "/components_reads";
        output_components_kmer_dir = output_path_str + "/components_kmer";
        comp_read_prefix = "/comp";
        comp_kmer_prefix = "/kmer";
        
        has_pair = false;
        
        duplicate_removed_read_dir = output_path_str + "/duplicate_remove_reads";
        
        boost::filesystem::path output_path(output_path_str);
        if( !boost::filesystem::exists(output_path_str) )
            add_directory(output_path);        
    }
    
    //for double stranded read
    Local_files(s &base_path_arg, s &output_dir_name,s &input_kmer_name, s &input_kmer_name_2,
                s &input_read_name, s &input_read_name_2, s &log_name, 
               s &contig_name_arg, s &comp_name_arg, s &kmer_name_arg) : 
               base_path_str(base_path_arg) , contig_name(contig_name_arg), 
               comp_name(comp_name_arg), kmer_name(kmer_name_arg)
    {
        output_path_str = base_path_str + output_dir_name;
        
        input_kmer_path = base_path_str + "/test_data" + input_kmer_name;
        input_kmer_path_2 = base_path_str + "/test_data" + input_kmer_name_2;
        
        input_read_path = base_path_str + "/test_data" + input_read_name; 
        input_read_path_2 = base_path_str + "/test_data" + input_read_name_2; 
        
        log_filename_path = base_path_str + "/log" + log_name;
        output_contig_path = output_path_str + contig_name;
        output_comp_path = output_path_str + comp_name;
        output_kmer_path = output_path_str + kmer_name;
        deleted_contig_path = output_path_str + "/deleted_kmer";  
        
        output_components_read_dir = output_path_str + "/components_reads";
        output_components_kmer_dir = output_path_str + "/components_kmer";
        comp_read_prefix = "/comp";
        comp_kmer_prefix = "/kmer";
        
        has_pair = true;
        
        duplicate_removed_read_dir = output_path_str + "/duplicate_remove_reads";
        
        boost::filesystem::path output_path(output_path_str);
        if( !boost::filesystem::exists(output_path_str) )
            add_directory(output_path);
    }
    
    Local_files& operator = (const Local_files & lf)
    {
        base_path_str = lf.base_path_str;
        input_kmer_path = lf.input_kmer_path;
        input_read_path = lf.input_read_path;    
        output_contig_path = lf.output_contig_path;        
        output_kmer_path = lf.output_kmer_path;
        deleted_contig_path = lf.deleted_contig_path;

        //common
        log_filename_path = lf.log_filename_path;
        output_comp_path = lf.output_comp_path;
        output_path_str = lf.output_path_str;

        output_components_read_dir = lf.output_components_read_dir;
        output_components_kmer_dir = lf.output_components_kmer_dir;
        comp_read_prefix = lf.comp_read_prefix;
        comp_kmer_prefix = lf.comp_kmer_prefix;
        
        //in case paired read
        input_kmer_path_2 = lf.input_kmer_path_2;
        input_read_path_2 = lf.input_read_path_2;  
        has_pair = lf.has_pair;
        
        duplicate_removed_read_dir = lf.duplicate_removed_read_dir;
    
        //input name
        contig_name = lf.contig_name;
        comp_name = lf.comp_name;
        kmer_name = lf.kmer_name;
    }
    
    void set_input_pair_path(s &input_kmer_name , s &input_kmer_name_2, 
                             s &input_read_name,  s &input_read_name_2)
    {
        input_kmer_path = base_path_str + "/test_data" + input_kmer_name;        
        input_kmer_path_2 = base_path_str + "/test_data" + input_kmer_name_2;
        input_read_path = base_path_str + "/test_data" + input_read_name; 
        input_read_path_2 = base_path_str + "/test_data" + input_read_name_2; 
        has_pair = true;
    }
    void set_if_use_pair(bool is_use_pair)
    {
        has_pair = is_use_pair;
    }
    void set_output_dir(s &output_dir_name)
    {
        output_path_str = base_path_str + output_dir_name;
        output_contig_path = output_path_str + contig_name;
        output_comp_path = output_path_str + comp_name;
        output_kmer_path = output_path_str + kmer_name;
        deleted_contig_path = output_path_str + "/deleted_kmer";  
        
        output_components_read_dir = output_path_str + "/components_reads";
        output_components_kmer_dir = output_path_str + "/components_kmer";
    }
    
    //for single
    std::string base_path_str;
    std::string input_kmer_path;
    std::string input_read_path;    
    std::string output_contig_path;        
    std::string output_kmer_path;
    std::string deleted_contig_path;
    
    std::string duplicate_removed_read_dir;
    
    // after contig graph dump
    std::string output_components_read_dir;
    std::string output_components_kmer_dir;        
    std::string comp_read_prefix;
    std::string comp_kmer_prefix;
    
    //common
    std::string log_filename_path;
    std::string output_comp_path;
    std::string output_path_str;
    
    //in case paired read
    std::string input_kmer_path_2;
    std::string input_read_path_2;  
    bool has_pair;
private:
    //input name
    std::string contig_name;
    std::string comp_name;
    std::string kmer_name;    
};

#endif	/* LOCAL_FILE_STRUCTURE_H */

