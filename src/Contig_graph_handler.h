/*
 * File:   Contig_graph_handler.h
 * Author: bx
 *
 * Created on September 2, 2017, 4:33 PM
 */

#ifndef CONTIG_GRAPH_HANDLER_H
#define	CONTIG_GRAPH_HANDLER_H

#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <stack>
#include <set>
#include <stdlib.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include "local_file_structure.h"

#include "shc_type.h"

#include "encoding.h"

#include "Contig_handler.h"
#include "Kmer_handler.h"
#include "shc_google_sparsehash.h"
#include <metis.h>
#include "json_parser.h"
#include "File_dumper.h"

#define ONLY_EDGE_WEIGHT "001"
#define METIS_NODE_OFFSET 1
#define EMPTY_KEY (std::pow(2,sizeof(contig_num_t)*8-1)-1)
#define DELETED_KEY (std::pow(2,sizeof(contig_num_t)*8-1)-1000)
#define METIS_IMBALANCE_PRECISION 1000.0
#define SIMPLE_COMPONENT "s_"
#define COMPLEX_COMPONENT "c_"
#define RE_PARTITION "r"

#define ANY_LETTER 'N'
#define ONE_MB (1024*1024)
#define TEN_MB (1024*1024*1024)

#define READ_PROGRESS_STEP 1000000
#define KMER_PROGRESS_STEP 1000000

extern uint8_t get_num_bit[256];

//#define LOG_CONTIG_GRAPH
//#define LOG_METIS

//#define ASSIGN_BEST // comment it if want assign all

struct bundled_contig_index {
    bundled_contig_index(const contig_num_t& p=IMPOSSIBLE_CONTIG_NUM):
                        contig_index(p) {}
    contig_num_t contig_index;
};

struct bundled_weight {
    bundled_weight(const contig_num_t& p=IMPOSSIBLE_CONTIG_NUM):
                        weight(p) {}
    contig_edge_weight_t weight;
};

class Contig_graph_handler
{
    friend class Contig_handler;
    friend class Kmer_handler;

    //IMPORTANT MUST USE VEC, since we implicitly assume vd are 0,1.. in metis
    typedef boost::adjacency_list < boost::vecS, boost::vecS,
            boost::undirectedS, bundled_contig_index, bundled_weight > graph_t;

    typedef boost::graph_traits < graph_t >::vertex_descriptor vd_t;
    typedef boost::graph_traits < graph_t >::vertex_iterator vi_t;

    typedef boost::graph_traits < graph_t >::edge_descriptor ed_t;
    typedef boost::graph_traits < graph_t >::edge_iterator ei_t;

    typedef std::pair<ei_t, ei_t> eip_t;

    //typedef tsl::hopscotch_map<contig_num_t, vd_t> Contig_vertex_map;
    //typedef tsl::hopscotch_map<contig_num_t, vd_t>::iterator Contig_vertex_map_iterator;
    //typedef tsl::hopscotch_map<contig_num_t, vd_t>::const_iterator Contig_vertex_map_const_iterator;

    struct Metis_input {
        Metis_input() {xadj.push_back(0);}
        // m : number of edges, n number of nodes, see metis manual page 23
        void metis_reset_data(int m, int n)
        {
            xadj.resize(n+1);
            xadj[0] = 0;
            adjncy.resize(2*m);
            adjwgt.resize(2*m);
        }
        void metis_deallocate()
        {
            std::vector<idx_t>().swap(adjncy);
            std::vector<idx_t>().swap(xadj);
            std::vector<idx_t>().swap(adjwgt);
        }
        std::vector<idx_t> adjncy;
        std::vector<idx_t> xadj;
        std::vector<idx_t> adjwgt;
    };

    struct Explorable_contig_set {
        Explorable_contig_set() {}
        Explorable_contig_set(contig_num_t num_contig_) :
            total_num_contig(num_contig_), explorable_contig(num_contig_, true)
        {
            std::cout << "num_contig " << num_contig_ << std::endl;
            shc_log_info(shc_logname, "num_contig %u\n", total_num_contig);
            explore_start = 0;
            num_explored_contig = 0;
        }

        inline contig_num_t get_set_next_explorable_contig()
        {
            contig_num_t i = total_num_contig+1;
            for(i=explore_start; i<total_num_contig; i++)
            {
                if(explorable_contig[i])
                {
                    explorable_contig[i] = false;
                    num_explored_contig++;
                    explore_start = i+1;
                    break;
                }
            }
            return (i);
        }

        inline void set_explored(contig_num_t p)
        {
            explorable_contig[p] = false;
            num_explored_contig++;
        }

        inline void erase_explored(contig_num_t p)
        {
            //shc_log_warning("before\n");
            explorable_contig[p] = true;
            //shc_log_warning("after\n");
            num_explored_contig--;

            if(p < explore_start)
                explore_start = p;
        }

        inline bool is_explorable(contig_num_t i) {return explorable_contig[i]; }
        inline contig_num_t get_num_explored_contig() {return num_explored_contig; }
        contig_num_t get_num_explorable_contig()
        {
            int diff= static_cast<int>(total_num_contig) -
                  static_cast<int>(num_explored_contig);
            if(diff < 0 )
            {
                shc_log_error("diff < 0\n");
                exit(1);
            }
            return total_num_contig - num_explored_contig;
        }

        bool is_explore_all() {return total_num_contig==num_explored_contig; }

        contig_num_t explore_start; //smallest_explorable_contig
        std::vector<bool> explorable_contig;
        contig_num_t total_num_contig;
        contig_num_t num_explored_contig;
    };

    struct Contig_vertex_map {
        Contig_vertex_map(contig_num_t num_contig_) : vd_array(num_contig_, num_contig_+1)
        {impossible_contig = num_contig_+1;}

        inline void set_contig_vd(contig_num_t i, vd_t vd)
        {
            assert(vd_array[i] == impossible_contig);
            vd_array[i] = vd;
        }

        // return true if there is element
        inline bool get_contig_vd(contig_num_t i, vd_t & vd)
        {
            if(vd_array[i]!=impossible_contig)
            {
                vd = vd_array[i];
                return true;
            }
            else
            {
                return false;
            }

        }
        contig_num_t impossible_contig;
        contig_num_t num_contig;
        std::vector<vd_t> vd_array;
    };

    struct Conn_contig_set {
        Conn_contig_set() {}
        Conn_contig_set(contig_num_t total_contig_) :
            num_contig(0), connectable_contig(total_contig_, true) {}
        inline void reset()
        {
            for(uint64_t i=0; i<num_contig; i++)
            {
                connectable_contig[connected_contig[i]] = true;
            }
            connected_contig.clear();
            num_contig = 0;
        }
        inline void set_conn_contig(contig_num_t i)
        {
            if(connectable_contig[i])
            {
                connectable_contig[i] = false;
                connected_contig.push_back(i);
                num_contig++;
            }
        }
        inline bool is_contig_connectable(contig_num_t i)
        {
            return connectable_contig[i];
        }
        inline bool is_contig_connected(contig_num_t i)
        {
            return !connectable_contig[i];
        }
        std::vector<bool> connectable_contig;
        std::vector<contig_num_t> connected_contig;
        contig_num_t num_contig;
    };


public:
    Contig_graph_handler(kmer_len_t length, Kmer_handler * khp,
         Contig_handler * chp, Shannon_C_setting & setting_):
         k1mer_len(length), kh(khp), ch(chp),
         component_array(chp->num_contig, IMPOSSIBLE_COMP_NUM),
         setting(setting_), num_test(setting_.contig_graph_setup.num_test),
         metis_setup(setting_.metis_setup), lf(setting_.local_files),
         kmer_length(khp->kmer_length), explorable_contig_set(chp->num_contig),
         is_double_stranded(setting_.is_double_stranded),
         fasta_dumper(INIT_SIZE, THRESH, setting_.has_single, setting_.has_pair),
         kmer_dumper(), conn_contig_set(chp->num_contig),
         is_assign_best(setting_.contig_graph_setup.is_assign_best)
    {
        //explorable_contig_set.set_empty_key(chp->num_contig+2);
        //explorable_contig_set.set_deleted_key(chp->num_contig+3);
        //explorable_contig_set.resize(chp->num_contig*2);

        //contig_stack.reserve(chp->num_contig);
        if(setting_.has_single)
            num_single_read_file = get_num_seq(setting_.local_files.input_read_path);
        if(setting_.has_pair)
            num_pair_read_file = get_num_seq(setting_.local_files.input_read_path_1);

        curr_component_num = 0;
        is_set_collect_comp_num = false;
        accum_collect_contig_num = 0;
        if(metis_setup.is_multiple_partition)
        {
            shc_log_info(shc_logname, "multi-parition\n");
            component_array_aux.assign(chp->num_contig, IMPOSSIBLE_COMP_NUM);
        }

        num_complex_comp = 0;

        non_partition_size = metis_setup.partition_size * 1; // change to 1 if want to keep the size the same
        first_n_graph = 10;
        give_up_num = 20;
    };

    void run_contig_graph_handler();
    void run_contig_graph_partition();
    void remove_read_duplicate(std::string out_dir_path);
    void group_components();
    void assign_reads_to_components(std::string input_read_path);
    void assign_reads_to_components_mmap(std::string input_read_path);
    void assign_paired_read_to_components(
                        std::string & read_1_path, std::string & read_2_path);
    void assign_paired_read_to_components_mmap(
                        std::string & read_1_path, std::string & read_2_path);
    void assign_kmer_to_components();
    void dump_component_array(std::string & filename);
    void dump_graph_to_metis_file();

    void load_data_and_partition_reads();

    int get_num_components() {return curr_component_num;}


    void log_metis_input_data();
    void log_contig_graph_edge(bool has_detail);
    void log_component_type();

private:
    void break_and_keep_component(graph_t & graph);
    bool assign_comp_without_graph();
    void assign_comp_with_graph();
    void get_connected_contig_index_without_graph( Kmer_counter_map_iterator & it,
          std::vector<contig_num_t> & contig_stack,
          char left_letter, char right_letter);
    void run_metis_and_assign_comp(graph_t & graph, std::vector<comp_num_t> & array,
                                    std::vector<idx_t> & part);
    void add_weight_to_cut_and_update_graph(graph_t & graph, std::vector<idx_t> * part);

    void get_connected_contig_index(Kmer_counter_map_iterator & it, graph_t & graph,
        Contig_vertex_map & contig_vertex_map, std::vector<contig_num_t> & contig_stack,
        vd_t curr_vd, char before_letter, char next_letter);
    void update_graph(Contig_vertex_map & contig_vertex_map, graph_t & graph, vd_t vd_i, contig_num_t j);

    void create_metis_array_input(graph_t & graph);

    void create_metis_file_format_from_graph(std::string & filename);

    void update_comp_set(std::set<comp_num_t> & comp_set, bool is_forward,
         std::vector<comp_num_t> & comp_array,  char * seq,  int len);

    void update_comp_map(std::map<comp_num_t, int> & comp_map, bool is_forward,
         std::vector<comp_num_t> & comp_array,  char * seq,  int len);

    void update_comp_count(std::vector<comp_num_t> & components,
                         std::vector<comp_num_t> & counts, bool is_forward,
         std::vector<comp_num_t> & comp_array,  char * seq,  int len);

    bool decide_best_comp(std::map<comp_num_t, int> & comp_count,
                                                        comp_num_t & best_comp);

    inline comp_num_t assign_single_best_comp(struct Single_dumper & mfs,
        std::map<comp_num_t, int> & comp_count, char * seq_ptr, char * header_ptr,
                                            int seq_len, bool is_forward);

    inline void assign_single_all_comp(struct Single_dumper & mfs,
        std::map<comp_num_t, int> & comp_count, char * seq_ptr, char * header_ptr,
                                            int seq_len, bool is_forward);

    inline comp_num_t assign_pair_best_comp(struct Single_dumper & mfs_p1,
        struct Single_dumper & mfs_p2, std::map<comp_num_t, int> & comp_count,
        char * seq_ptr_p1, char * seq_ptr_p2, char * header_ptr_p1,
        char * header_ptr_p2, int seq_len_p1, int seq_len_p2, bool is_forward);
    inline void assign_pair_all_comp(struct Single_dumper & mfs_p1,
        struct Single_dumper & mfs_p2, std::map<comp_num_t, int> & comp_count,
        char * seq_ptr_p1, char * seq_ptr_p2, char * header_ptr_p1,
        char * header_ptr_p2, int seq_len_p1, int seq_len_p2, bool is_forward);

    inline void assign_single_read_to_file_mmap_helper(struct Single_dumper & mfs,
              char * header_ptr, char * seq_ptr, int seq_len, bool is_forward);
    inline void assign_pair_read_to_file_mmap_helper(struct Single_dumper & mfs_p1,
               struct Single_dumper & mfs_p2, char * header_ptr_p1,
               char * header_ptr_p2, char * seq_ptr_p1, char * seq_ptr_p2,
               int seq_len_p1, int seq_len_p2, bool is_forward);

    void assign_single_read_to_file(std::string & header, std::string & sequence,
            std::vector<std::shared_ptr<std::ofstream> > & files);
    void assign_single_read_to_file_mmap(struct Single_dumper & mfs,
                                   char * header_ptr, char * seq_ptr, int len);

    void assign_pair_read_to_file_mmap(struct Single_dumper & mmap_files_p1,
                                       struct Single_dumper & mmap_files_p2,
                                   char * header_ptr_p1, char * header_ptr_p2,
                                   char * seq_ptr_p1,    char * seq_ptr_p2,
                                   int len_p1,           int len_p2);

    void assign_pair_read_to_file(std::string & header1, std::string & sequence1,
                   std::string & header2, std::string & sequence2,
             std::vector<std::shared_ptr<std::ofstream> > & files1,
                          std::vector<std::shared_ptr<std::ofstream> > & files2);

    void count_component_size();
    void log_comp_content();
    void log_contigs(std::string file_path);

    kmer_len_t k1mer_len;

    Explorable_contig_set explorable_contig_set;
    std::vector<contig_num_t> contig_stack;
    Conn_contig_set conn_contig_set;

    struct FASTA_dumper fasta_dumper;
    struct KMER_dumper kmer_dumper;


    // capture mapping from contig index to component index, the comp in
    // component_array and its aux are interleaved. i.e, no first k comp in one
    // vector and last k in auxiliary vector.
    std::vector<comp_num_t> component_array; // index represents contig id, content stores comp id
    std::vector<comp_num_t> component_array_aux; // used for re-partition

    std::vector<std::string> component_type;
    int num_complex_comp;

    comp_num_t curr_component_num; // after processing the largest comp num is this number -1

    comp_num_t collect_comp_num;
    int accum_collect_contig_num;
    bool is_set_collect_comp_num;

    Kmer_handler * kh;
    Contig_handler * ch;

    struct Metis_input metis_input;

    struct Shannon_C_setting & setting;
    struct Metis_setup metis_setup;
    Local_files & lf;
    uint8_t kmer_length;
    int num_test;
    bool is_double_stranded;
    int num_read_files;
    int num_kmer_files;

    uint64_t num_single_read_file;
    uint64_t num_pair_read_file;

    uint64_t total_query;
    std::vector<uint64_t> component_size;

    bool is_assign_best;

    idx_t non_partition_size;
    uint64_t first_n_graph;
    int give_up_num;
};

#endif	/* CONTIG_GRAPH_HANDLER_H */
