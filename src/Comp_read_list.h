/*
 * File:   Comp_read_list.h
 * Author: bx
 *
 * Created on October 7, 2017, 1:16 PM
 */

#ifndef COMP_READ_LIST_H
#define	COMP_READ_LIST_H

#include <vector>
#include <string.h>
#include <fstream>
#include "shc_type.h"
#include "local_file_structure.h"
#include "log.h"


struct Read_acc {
    Read_acc (char * s_p, read_length_t len_p):
            read_ptr(s_p), len(len_p) {}
    Read_acc & operator = (const Read_acc & l)
    {
        read_ptr = l.read_ptr;
        len = l.len;
    }
    char * read_ptr;
    read_length_t len;
};

struct Comp_read_list {
    Comp_read_list(bool is_compress_) :
        is_compress(is_compress_), read_length(0)
    {
        read_len_map.set_empty_key(IMPOSSIBLE_READ_NUM-1);
        read_len_map.set_deleted_key(IMPOSSIBLE_READ_NUM-2);
        num_read = 0;
    }

    void setup(const std::string read_path)
    {
        get_read_length(read_path);
        size_t estimated_num_read = estimate_num_read(read_path, read_length);
#ifdef SHOW_SEQ_GRAPH_PROGRESS
        std::cout << "estimated_num_read is " << estimated_num_read << std::endl;
#endif
        read_list.reserve(estimated_num_read*read_length);
    }

    void add_read(std::string & base)
    {
        //if this read is truncated
        assert(base.size() <= read_length);
        if(base.size() < read_length)
        {
            read_len_map[num_read] = base.size();
            base.insert(base.end(), read_length-base.size(), 'C');
        }

        if(!is_compress)
        {
            read_list.insert(read_list.end(), base.begin(), base.end());
        }
        else
        {
            encode_in_place((char*)&base.at(0),(uint8_t*)&base.at(0), read_length);
            read_list.insert(read_list.end(), &base.at(0),
                                                 &base.at(compress_length));
        }
        num_read++;
    }

    Read_acc get_read(read_num_t i)
    {
        if(!is_compress)
        {
            Read_len_map_iterator it = read_len_map.find(i);
            if( it == read_len_map.end())
                return Read_acc((char*)&(read_list.at(i*read_length)), read_length);
            else
                return Read_acc((char*)&(read_list.at(i*read_length)), it->second);

        }
        else
        {
            if(curr_read.size()==0)
                curr_read.resize(read_length);

            Read_len_map_iterator it = read_len_map.find(i);
            if( it == read_len_map.end())
            {
                if(i != curr_read_i)
                {
                    curr_read_i = i;
                    decode_byte_list(&read_list.at(i*compress_length),
                            (char *)&curr_read.at(0), read_length);
                }
                return Read_acc((char*)&curr_read.at(0), read_length);
            }
            else
            {
                if(i != curr_read_i)
                {
                    curr_read_i = i;
                    decode_byte_list(&read_list.at(i*compress_length),
                                (char *)&curr_read.at(0), it->second);
                }
                return Read_acc((char*)&curr_read.at(0), it->second);
            }
        }
    }

    void get_read_length(const std::string & filename)
    {
        std::string trunc_symbol("trunc");
        std::string header, read_base;
        std::ifstream file_reader(filename.c_str());
        while (std::getline(file_reader, header) &&
            std::getline(file_reader, read_base) )
        {
            //if no match is found
            if(header.find(trunc_symbol) == std::string::npos)
            {
                read_length = read_base.size();
                compress_length = static_cast<read_length_t>
                        (std::ceil(static_cast<double>(read_length)/4.0));
                break;
            }
        }
        std::cout << "read length is " << read_length <<std::endl;
        if(header.empty() || read_base.empty())
        {
            shc_log_error("%s file contain empty contents\n", filename.c_str());
            exit(1);
        }
    }

    void deallocate_resource()
    {
        std::vector<uint8_t>().swap(read_list);
        std::string().swap(curr_read);
    }

    read_num_t get_num_read(){return num_read;}
    read_num_t get_read_id(){return num_read;}
    read_length_t get_read_length() {return read_length;}

private:
    read_num_t num_read;
    read_length_t read_length;
    read_length_t compress_length;
    bool is_compress;

    std::vector<uint8_t> read_list;
    std::vector<uint8_t> read_mate_type;
    std::vector<read_num_t> mate_read_num;

    Read_len_map read_len_map;

    std::string curr_read;
    read_num_t curr_read_i;  // allow caching
};


#endif	/* COMP_READ_LIST_H */
