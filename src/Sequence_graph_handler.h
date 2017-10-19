/*
 * File:   Sequence_graph_handler.h
 * Author: bx
 *
 * Created on October 2, 2017, 8:15 PM
 */

#ifndef SEQUENCE_GRAPH_HANDLER_H
#define	SEQUENCE_GRAPH_HANDLER_H


#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include <vector>
#include <string.h>
#include <iostream>
#include <fstream>
#include <set>
#include <deque>

#include "shc_google_sparsehash.h"
#include "json_parser.h"
#include "Kmer_handler.h"
#include "local_file_structure.h"
#include "log.h"
#include "Comp_read_list.h"
#include <cmath>

#define EMPTY_READ_NUM (std::pow(2, sizeof(read_num_t)*8)-1)
#define EMPTY_READ_LENGTH (std::pow(2, sizeof(read_length_t)*8)-1)
#define TWO_POWER_16 65536
#define TWO_POWER_8 256

class Sequence_graph_handler {
    friend class Kmer_handler;
    friend class Multibridge_handler;

    struct Bdg_read_info {
        Bdg_read_info() : read_id(IMPOSSIBLE_READ_NUM-10),
                                        start(IMPOSSIBLE_READ_LENGTH-10){}
        Bdg_read_info(read_num_t read_id_, read_length_t start_) :
                                        read_id(read_id_), start(start_){}

        Bdg_read_info (const Bdg_read_info & b)
        {
            read_id = b.read_id;
            start = b.start;
        }
        Bdg_read_info & operator = (const Bdg_read_info & b)
        {
            read_id = b.read_id;
            start = b.start;
            return (*this);
        }
        read_num_t read_id;
        read_length_t start;
    };

    struct bundled_node_p {
        bundled_node_p(const std::string & seq_="", const kmer_count_t &
                  count_=IMPOSSIBLE_KMER_NUM) : seq(seq_), count(count_) {}
        read_length_t seq_len(){return seq.size();}
        std::string seq;
        kmer_count_t count;
        std::vector<Bdg_read_info> reads_info;
    };

    struct bundled_edge_p {
        bundled_edge_p(const edge_weight_t & weight_=IMPOSSIBLE_EDGE_WEIGHT_NUM,
                       const kmer_count_t & copy_count_=0)
                            : weight(weight_), copy_count(copy_count_) {}
        edge_weight_t weight;
        kmer_count_t copy_count;
    };

public:
    typedef boost::adjacency_list<boost::listS, boost::listS,
            boost::bidirectionalS, bundled_node_p, bundled_edge_p> seq_graph_t;

    typedef boost::graph_traits < seq_graph_t >::vertex_descriptor vd_t;
    typedef boost::graph_traits < seq_graph_t >::vertex_iterator vi_t;

    typedef boost::graph_traits < seq_graph_t >::edge_descriptor ed_t;
    typedef boost::graph_traits < seq_graph_t >::edge_iterator ei_t;

    typedef boost::graph_traits<seq_graph_t>::in_edge_iterator in_ei_t;
    typedef boost::graph_traits<seq_graph_t>::out_edge_iterator out_ei_t;

    typedef std::pair<in_ei_t, in_ei_t> in_eip_t;
    typedef std::pair<out_ei_t, out_ei_t> out_eip_t;
    typedef std::pair<ei_t, ei_t> eip_t;

    typedef std::pair<vi_t, vi_t> vip_t;
    typedef std::pair<vd_t, vd_t> vdp_t;

    typedef std::pair<ed_t, bool> aer_t; // an edge reference

    //Google sparsehash
    typedef google::dense_hash_map<std::string, vd_t,
                                std::hash<std::string>, eqstr> Kmer_Node_map;
    typedef google::dense_hash_map<std::string, vd_t,
             std::hash<std::string>, eqstr>::iterator Kmer_Node_map_iterator;

    typedef google::dense_hash_map<read_num_t, std::vector<vd_t>,
                                    hash_u32, equ32> Read_path_map;
    typedef google::dense_hash_map<read_num_t, std::vector<vd_t>,
                        hash_u32, equ32>::iterator Read_path_map_iterator;

    typedef std::set<vd_t> Node_set_t;
    typedef std::set<vd_t>::iterator Node_set_iterator;

    struct Node_index {
        Node_index(vd_t vd_, read_length_t start_) : vd(vd_), start(start_) {}
        vd_t vd;
        read_length_t start;
    };

    Sequence_graph_handler(Kmer_handler * kh_, Shannon_C_setting & setting_,
                           uint64_t comp_i, bool is_compress, int max_hop_) :
                           kh(kh_), setting(setting_), read_list(is_compress),
                           max_hop(max_hop_), component_i(comp_i)
    {
        kmer_path = setting.local_files.output_components_kmer_dir
                + setting.local_files.comp_kmer_prefix + std::to_string(comp_i);
        read_path = setting.local_files.output_components_read_dir
                + setting.local_files.comp_read_prefix + std::to_string(comp_i);

        kmer_length = kh->kmer_length;

        read_list.setup(read_path);
    }

    Sequence_graph_handler(uint8_t kmer_length_, Shannon_C_setting & setting_,
                           uint64_t comp_i, bool is_compress, int max_hop_) :
                          kmer_length(kmer_length_), setting(setting_),
                          read_list(is_compress), max_hop(max_hop_),
                          component_i(comp_i)
    {
        kmer_path = setting.local_files.output_components_kmer_dir
                + setting.local_files.comp_kmer_prefix + std::to_string(comp_i);
        read_path = setting.local_files.output_components_read_dir
                + setting.local_files.comp_read_prefix + std::to_string(comp_i);

        read_list.setup(read_path);

        kmer_node_map.set_empty_key("");
        read_known_path_map.set_empty_key(IMPOSSIBLE_READ_NUM);
    }

    void * run_seq_graph_handler_helper();

    static void * run_seq_graph_handler(void * context)
    {
        return ((Sequence_graph_handler *)context)->run_seq_graph_handler_helper();
    }

    void load_kmer_and_build_graph();
    void build_kmer_graph_from_reads();

    void assign_read_to_xnode();
    bool is_read_bridge_node(std::string & read_seq, read_length_t start, vd_t vd);
    bool is_read_bridge_node(const char * read_ptr, read_length_t len,
                                        read_length_t info_start, vd_t vd);
    void node_link_read(vd_t vd, read_num_t read_id, read_length_t start);

    void condense_graph();

    vd_t merge_two_vertices(vd_t vd1, vd_t vd_2);

    void update_xnode_set();
    bool is_xnode(vd_t vd);

    void bridge_all_xnodes();
    void refresh_bridging_reads(vd_t vd);
    bool is_read_bridge_left(Bdg_read_info & read_info, vd_t vd);
    bool is_read_bridge_right(Bdg_read_info & read_info, vd_t vd);
    bool is_xnode_bridged(vd_t vd);

    void perform_bridging(vd_t vd);
    bool is_read_bridge_node(read_num_t read_id, read_length_t start, vd_t vd);
    vd_t create_left_node(ed_t ed);
    vd_t create_right_node(ed_t ed);
    void local_condense_bridged(std::vector<vd_t> & u_list, std::vector<vd_t> & w_list);
    void local_condense(vd_t vd, std::set<vd_t> & vd_remove_set);
    void condense_left(vd_t vd, vd_t vd_new, std::set<vd_t> & vd_remove_set);
    void condense_right(vd_t vd, vd_t vd_new, std::set<vd_t> & vd_remove_set);
    bool condense(vd_t vd, bool us_check_left , vd_t vd_new);
    bool is_condensable(vd_t vd, bool is_left);

    void break_self_loops();

    void break_all_cycles();
    bool find_cycle(std::set<vd_t> & acyclic_node_set, std::deque<vd_t> & cycle_path);
    bool is_node_inside_cycle(vd_t vd, std::set<vd_t> & acyclic_node_set, std::deque<vd_t> & cycle_path);
    void break_cycle(std::deque<vd_t> & cycle_path);
    void break_cycle_splitting(vd_t vd_split, vd_t vd_c_in, vd_t vd_c_out,
        int sum_copy_count, int sum_c_in_and_o_out);


    void find_known_path(int curr_hop);
    bool is_partial_match(const char * ptr1, const read_length_t len1,
                          const char * ptr2, const read_length_t len2);
    bool search_seq(const char * seq_ptr, read_length_t seq_len,
                    vd_t vd, read_length_t start,
                    int curr_hop, std::vector<vd_t> & path);

    void find_mate_pairs();

    void log_all_node(bool is_log_node, bool is_log_read);
    void log_xnode_set(bool is_log_node, bool is_log_read);
    void log_read_seq(Read_acc acc);
    void log_node_info(vd_t vd, bool is_log_nodes , bool is_log_read_info);
    void log_edge_weight();
    void log_seq_only(vd_t vd, read_length_t start);


private:
    //for each component, needs following data
    seq_graph_t graph;
    Comp_read_list read_list;

    Kmer_Node_map kmer_node_map;
    Node_set_t xnode_set;

    //std::set<std::vector<vd_t> > known_path_set;
    Read_path_map read_known_path_map;
    std::set<vdp_t> mate_nodes_set;

    uint8_t kmer_length;
    Kmer_handler * kh;
    Shannon_C_setting & setting;
    int max_hop;

    std::string kmer_path;
    std::string read_path;

    Block_timer timer;
    uint64_t component_i;
};
#endif	/* SEQUENCE_GRAPH_HANDLER_H */
