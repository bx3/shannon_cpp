/*
 * File:   Comp_read_list.h
 * Author: bx
 *
 * Created on October 7, 2017, 1:16 PM
 */

#ifndef SINGLE_READ_LIST_H
#define	SINGLE_READ_LIST_H

#include <vector>
#include <string.h>
#include <fstream>
#include "shc_type.h"
#include "local_file_structure.h"
#include "log.h"
#include "json_parser.h"

#define SINGLE 0
#define PAIR_1 1
#define PAIR_2 2
#define SUPERR 3
#define OVERLAP_CHAR 'L'


struct Read_acc {
    Read_acc (char * s_p, read_length_t len_p, read_count_t read_count_p,
                            double read_prob_):
            read_ptr(s_p), len(len_p), read_count(read_count_p),
            read_prob(read_prob_) {}
    Read_acc (const Read_acc & a)
    {
        read_ptr = a.read_ptr;
        len = a.len;
        read_count = a.read_count;
        read_prob = a.read_prob;
    }
    Read_acc & operator = (const Read_acc & l)
    {
        read_ptr = l.read_ptr;
        len = l.len;
        read_count = l.read_count;
        read_prob = l.read_prob;
        return *this;
    }

    inline read_count_t get_estimated_read_count() {return std::ceil(read_count);} // without prob /read_prob

    char * read_ptr;
    read_length_t len;
    read_count_t read_count;
    double read_prob;
};

struct Read_acc_pair {
    Read_acc_pair(const Read_acc & a1, const Read_acc & a2):
            acc1(a1), acc2(a2) {}
    struct Read_acc acc1;
    struct Read_acc acc2;
};

struct Single_read_list {
    typedef google::dense_hash_map<std::string, read_num_t,
                   std::hash<std::string>, eqstr> Read_Index_map;
    typedef google::dense_hash_map<std::string, read_num_t,
                std::hash<std::string>, eqstr>::iterator Read_Index_map_iterator;

    Single_read_list(read_length_t read_length_, bool is_compress_,
                    Shannon_C_setting & setting_, uint8_t read_type_):
        read_length(read_length_), is_compress(is_compress_),
        setting(setting_), read_type(read_type_)
    {
        compress_length = static_cast<read_length_t>
                            (std::ceil(static_cast<double>(read_length)/4.0));
        total_num_read = 0;
        num_read = 0;
        read_index_map.set_deleted_key("d");
        read_index_map.set_empty_key("");
        //curr_read_i = IMPOSSIBLE_READ_NUM;
    }

    void setup(size_t total_num_read_)
    {
#ifdef SHOW_SEQ_GRAPH_PROGRESS
        std::cout << "estimated_num_read is " << estimated_num_read << std::endl;
#endif
        //total_num_read = total_num_read_;
        read_list.reserve(total_num_read_*read_length);
        read_len_map.assign(total_num_read_, 0);
        if(read_type != SINGLE)
            read_pair_index.reserve(total_num_read_);
    }

    // return if it is a new read, return on argument the local read index
    bool add_read(std::string & base, read_num_t & read_index, double read_prob)
    {
        total_num_read++;

        if(base.size() > read_length)
        {
            base.resize(read_length);
        }

        Read_Index_map_iterator it = read_index_map.find(base);
        if(it == read_index_map.end())
        {
            read_index = num_read;
            read_index_map.insert(std::make_pair(base, read_index));
            local_read_index.push_back(read_index);
            read_counts.push_back(1);
            read_probs.push_back(read_prob);

            if(base.size() < read_length)
            {
                read_len_map[read_index] = base.size();
                base.insert(base.end(), read_length-base.size(), 'T');
            }
            num_read++;

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
            return true;
        }
        else
        {
            read_index = it->second;
            read_counts[read_index]++;
            //read_probs[read_index] += read_prob; since if read are equal prob are equal
            local_read_index.push_back(read_index);
            return false;
        }

    }

    // assume pair read are processed at the same time
    void link_read(read_num_t other_read_index)
    {
        read_pair_index.push_back(other_read_index);
    }

    read_num_t get_pair_read_index(read_num_t my_read_index)
    {
        return read_pair_index[my_read_index];
    }

    read_count_t get_count(read_num_t i)
    {
        if (i>= num_read)
        {
            shc_log_error("comp %d, Access impossible read %u\n", comp_id, i);
            exit(1);
        }
        return read_counts[i];
    }

    void correct_read_count(std::string & base, read_count_t & read_count)
    {
        Read_Index_map_iterator it = read_index_map.find(base);
        assert(it != read_index_map.end());
        read_counts[it->second] = read_count;
    }

    void correct_read_count(read_num_t i, read_count_t read_count)
    {
        read_counts[i] = read_count;
    }

    // input read index local to this collection
    Read_acc get_read(read_num_t i)
    {
        if (i>= num_read)
        {
            shc_log_error("comp %d, Access impossible read %u\n", comp_id, i);
            exit(1);
        }

        //if(!is_compress)
        //{
            if(read_len_map[i]==0)
            {
                return Read_acc((char*)&(read_list.at(i*read_length)),
                            read_length, read_counts[i], read_probs[i]);
            }
            else
            {
                return Read_acc((char*)&(read_list.at(i*read_length)),
                            read_len_map[i], read_counts[i], read_probs[i]);
            }

            //Read_len_map_iterator it = read_len_map.find(i);
            //if( it == read_len_map.end())
            //    return Read_acc((char*)&(read_list.at(i*read_length)),
            //                read_length, read_counts[i]);
            //else
            //    return Read_acc((char*)&(read_list.at(i*read_length)),
            //                it->second, read_counts[i]);
        //}
        /*
        else
        {
            if(curr_read.size()==0)
                curr_read.resize(read_length);

            Read_len_map_iterator it = read_len_map.find(i);
            if( it == read_len_map.end())
            {
                //if(i != curr_read_i)
                //{
                //    curr_read_i = i;
                    decode_byte_list(&read_list.at(i*compress_length),
                            (char *)&curr_read.at(0), read_length);
                //}
                return Read_acc((char*)&curr_read.at(0), read_length, read_counts[i]);
            }
            else
            {
                //if(i != curr_read_i)
                //{
                //    curr_read_i = i;
                    decode_byte_list(&read_list.at(i*compress_length),
                                (char *)&curr_read.at(0), it->second);
                //}
                return Read_acc((char*)&curr_read.at(0), it->second, read_counts[i]);
            }
        }
        */
    }

    void deallocate_resource()
    {
        std::vector<uint8_t>().swap(read_list);
        std::string().swap(curr_read);
    }

    void deallocate_unneeded_resource_after_loading_read()
    {
        read_index_map.clear();
    }

    void clear()
    {
        total_num_read = 0;
        num_read = 0;
        read_list.clear();
        read_list.shrink_to_fit();
        read_len_map.clear();
        curr_read.clear();
        local_read_index.clear();
        curr_read_i = IMPOSSIBLE_READ_NUM;
        comp_id = -1;
    }
    //from align index to local index
    inline read_num_t get_local_read_index(read_num_t i){return local_read_index[i];}

    inline read_num_t get_total_num_read(){return total_num_read;}
    inline read_num_t get_num_reads() {return num_read;}
    inline read_num_t get_read_id(){return total_num_read;}
    inline read_length_t get_read_length() {return read_length;}

    inline void set_comp_id(int i) {comp_id = i;}

    inline read_num_t get_num_read() {return  num_read; }

private:
    uint8_t read_type;
    read_num_t total_num_read;  //or number of pairs
    read_num_t num_read;

    std::vector<uint8_t> read_list; //now only store unique read
    std::vector<read_num_t> read_pair_index;
    // store mapping from align read index to local read index
    std::vector<read_num_t> local_read_index;
    std::vector<read_count_t> read_counts;
    std::vector<double> read_probs;

    std::vector<read_length_t> read_len_map;
    std::string curr_read;
    read_num_t curr_read_i;  // allow caching

    //once set, should not be changed
    read_length_t read_length;
    bool is_compress;
    read_length_t compress_length;

    Read_Index_map read_index_map;

    Shannon_C_setting & setting;
    int comp_id;
};


#endif	/* COMP_READ_LIST_H */
