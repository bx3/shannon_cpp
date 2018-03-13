#ifndef COLLECT_READS
#define COLLECT_READS
#include "Single_read_list.h"

struct Super_read {
    Super_read(){}
    Super_read(read_num_t start_id, read_num_t stop_id, std::string &middle_):
            start_read_id(start_id), stop_read_id(stop_id), middle(middle_) {}
    Super_read & operator = (const Super_read & b)
    {
        start_read_id = b.start_read_id;
        stop_read_id = b.stop_read_id;
        middle = b.middle;
        return (*this);
    }
    bool operator == (const Super_read & sr) const
    {
        return start_read_id==sr.start_read_id && stop_read_id==sr.stop_read_id && middle==sr.middle;
    }
    read_num_t start_read_id;  // lcoal read index
    read_num_t stop_read_id;
    std::string middle;
};

struct Collect_reads {
    Collect_reads (bool has_pair_, bool has_single_,
                     read_length_t s_length_,
                     read_length_t p1_length_,
                     read_length_t p2_length_,
                     bool is_compress_,
                     Shannon_C_setting & setting_) :
                has_paired_read(has_pair_), has_single_read(has_single_),
                s_length(s_length_),
                s_reads(s_length_, is_compress_, setting_, SINGLE),
                p1_reads(p1_length_, is_compress_, setting_, PAIR_1),
                p2_reads(p2_length_, is_compress_, setting_, PAIR_2),
                setting(setting_)
    {
        total_num_read = 0;
        pair_read_start_num = 0;
        if(has_single_ && has_pair_)
        {
            L = (s_reads.get_read_length() + p1_reads.get_read_length() +
            p2_reads.get_read_length())/3;
        }
        else if (has_single_)
        {
            L = s_reads.get_read_length();
        }
        else
            L = (p1_reads.get_read_length() + p2_reads.get_read_length())/2;
    }

    //number of reads
    read_num_t get_num_single_read() {return s_reads.get_num_reads();}
    read_num_t get_num_paired_read()
    {
        assert(p1_reads.get_num_reads() == p2_reads.get_num_reads());
        return p1_reads.get_num_reads();
    }

    typedef std::map<read_num_t, Super_read> Read_sr_map;
    typedef std::map<read_num_t, Super_read>::iterator Read_sr_map_iterator;

    read_num_t get_total_num_reads()
    {

    }

    void clear()
    {
        s_reads.clear();
        p1_reads.clear();
        p2_reads.clear();
        total_num_read = 0;
        pair_read_start_num = 0;
        read_sr_map.clear();
        comp_id = -1;
        //std::string impossible_middle = "ZZ";
        //last_sr = Super_read(0, 0, impossible_middle);
    }

    //add reads
    //void add_single_read(std::string & read){ s_reads.add_read(read);  }
    //void add_paired_read(std::string & read1, std::string & read2)
    //                {p1_reads.add_read(read1); p2_reads.add_read(read2); }

    bool is_p1_read(read_num_t i) {return (i-s_reads.get_num_reads())%2 == 0;}
    bool is_sr(read_num_t virtual_index)
            { return read_sr_map.find(virtual_index)!= read_sr_map.end();}
    //there is a virtual mapping from j to reads, the order is
    // single reads , pair1 read pair2 read, interleaved

    //convert read index approriate to their local index
    //used for iteration
    Read_acc get_read(read_num_t i)
    {
        Read_sr_map_iterator it = read_sr_map.find(i);
        if(it != read_sr_map.end())
        {
            //if (i<s_reads.get_num_reads()+p1_reads.get_num_reads()) // ignore the p2 part
                return get_super_read(it->second);

        }

        if(i<s_reads.get_num_reads())
        {
            //shc_log_info(shc_logname, "%d read, single\n", i);
            return s_reads.get_read(i);
        }
        else if(i < vs_p2_read)
        {
            //shc_log_info(shc_logname, "%d read, pair 1\n", i);
            return p1_reads.get_read(i-vs_p1_read);
        }
        else if(i < get_num_reads())
        {
            //shc_log_info(shc_logname, "%d read, pair 2\n", i);
            return p2_reads.get_read(i-vs_p2_read);
        }
        else
        {
            shc_log_error("comp %d, encounter unknown num read %u\n", comp_id, i);
            exit(1);
        }
    }

    read_num_t get_num_reads()
    {
        size_t total_num =
               static_cast<size_t>(s_reads.get_num_reads()) +
               static_cast<size_t>(p1_reads.get_num_reads()) +
               static_cast<size_t>(p2_reads.get_num_reads()) ;
        if(total_num > IMPOSSIBLE_READ_NUM)
        {
            shc_log_error("comp %d, read num type unable to hold all read ids\n", comp_id);
            exit(1);
        }
        return total_num;
    }

    void link_two_reads(read_num_t read_index1, read_num_t read_index2)
    {
        p1_reads.link_read(read_index2);
        p2_reads.link_read(read_index1);
    }

    read_num_t convert_local_index_to_virtual_index(read_num_t i, bool is_p1)
    {
        if(is_p1)
        {
            return i + s_reads.get_num_reads();
        }
        else
        {
            return i + s_reads.get_num_reads() + p1_reads.get_num_reads();
        }
    }

    void log_all_super_read()
    {
        shc_log_info(shc_logname, "Start log all Super read\n");
        for(Read_sr_map_iterator it =read_sr_map.begin();
                                 it!=read_sr_map.end(); it++)
        {
            Read_acc acc = get_super_read(it->second);
            std::string sr(acc.read_ptr, acc.len);
            shc_log_info(shc_logname, "len %d: %s\n",sr.size(), sr.c_str());
        }
        shc_log_info(shc_logname, "Finish log all Super read\n");
    }

    void declare_read_finish()
    {
        vs_s_read = 0; //vs for virtual start
        vs_p1_read = s_reads.get_num_reads();
        vs_p2_read = s_reads.get_num_reads() + p1_reads.get_num_reads();
    }

    read_num_t convert_align_index_to_virtual_index(read_num_t i, bool is_p1)
    {
        if(is_p1)
        {
            return convert_local_index_to_virtual_index(
                                       p1_reads.get_local_read_index(i), true);
        }
        else
        {
            return convert_local_index_to_virtual_index(
                                       p2_reads.get_local_read_index(i), false);
        }
    }

    read_num_t convert_align_index_to_local_index(read_num_t i, bool is_p1)
    {
        if(is_p1)
            return p1_reads.get_local_read_index(i);
        else
            return p2_reads.get_local_read_index(i);
    }


    void add_super_read(read_num_t virtual_start_read_id, read_num_t virtual_stop_read_id,
            read_num_t l_start_read_id, read_num_t l_stop_read_id,
            std::string & middle)
    {
        Super_read sr(l_start_read_id, l_stop_read_id, middle);
        read_sr_map.insert(std::make_pair(virtual_start_read_id, sr));
        read_sr_map.insert(std::make_pair(virtual_stop_read_id, sr));
    }


    Read_acc get_super_read(Super_read & sr)
    {
        //shc_log_info(shc_logname, "get super read\n");
        // preparecurr_sr
        //if(sr == last_sr)
        //{
        //    return Read_acc(&curr_sr.at(0), curr_sr.size());
        //}
        //else
        //{
            Read_acc ra1 = p1_reads.get_read(sr.start_read_id);
            Read_acc ra2 = p2_reads.get_read(sr.stop_read_id);

            //shc_log_info(shc_logname, "ra1 %d, ra2 %d\n", ra1.len, ra2.len);

            curr_sr.resize(ra1.len);
            memcpy(&curr_sr.at(0), ra1.read_ptr, ra1.len);
            if(sr.middle.size() == 0)
            {
                curr_sr.resize(curr_sr.size()+ra2.len);
                memcpy(&curr_sr.at(ra1.len), ra2.read_ptr, ra2.len);
            }
            else if(sr.middle.at(0) != OVERLAP_CHAR)
            {
                curr_sr += sr.middle;
                curr_sr.resize(curr_sr.size()+ra2.len);
                memcpy(&curr_sr.at(ra1.len + sr.middle.size()), ra2.read_ptr, ra2.len);
            }
            else
            {
                //shc_log_info(shc_logname, "overlap super read\n");
                read_num_t overlap = std::stoi(sr.middle.substr(1));
                curr_sr.resize(curr_sr.size()+ra2.len - overlap);
                memcpy(&curr_sr.at(ra1.len), ra2.read_ptr+ overlap, ra2.len-overlap);
            }

            read_count_t sr_count = (ra1.read_count + ra2.read_count)/2;
            last_sr = sr;
            return Read_acc(&curr_sr.at(0), curr_sr.size(), sr_count);
        //}
    }

    void set_comp_id(int i)
    {
        comp_id = i;
        s_reads.set_comp_id(i);
        p1_reads.set_comp_id(i);
        p2_reads.set_comp_id(i);
    }
    read_num_t get_single_reads_num() {return s_reads.get_num_reads();}
    read_num_t get_p1_reads_start_num() {return s_reads.get_num_reads();}
    read_num_t get_p2_reads_start_num()
                {return s_reads.get_num_reads()+p1_reads.get_num_reads();}

    //data structures
    struct Single_read_list s_reads;
    struct Single_read_list p1_reads;
    struct Single_read_list p2_reads;

    read_num_t vs_s_read; //vs for virtual start
    read_num_t vs_p1_read;
    read_num_t vs_p2_read;

    read_num_t total_num_read;
    read_num_t pair_read_start_num;

    Read_sr_map read_sr_map; // sr for super read

    Super_read last_sr;
    std::string curr_sr;

    //setting parameter, not changed once set
    read_length_t s_length;
    bool has_paired_read;
    bool has_single_read;
    int comp_id;

    Shannon_C_setting & setting;
    read_length_t L;
};
#endif
