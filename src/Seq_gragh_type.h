/*
 * File:   Seq_gragh_type.h
 * Author: bx
 *
 * Created on November 14, 2017, 11:37 PM
 */
#ifndef SEQ_GRAGH_TYPE_H
#define	SEQ_GRAGH_TYPE_H

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/any.hpp>

#include <vector>
#include <string.h>

#include "shc_type.h"

#define NO_PAIR 0
#define PAIR_END_NODE 1
#define PAIR_FRONT_NODE 2

#define WHITE 0
#define GRAY  1
#define BLACK 2

typedef std::pair<char, char> Adj_bases; //left right base

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
    bool operator< (const struct Bdg_read_info & read_info) const
    {
        if(this->read_id < read_info.read_id)
            return true;
        else if (this->read_id > read_info.read_id)
            return false;
        else
        {
            if(this->start < read_info.start)
                return true;
            else
                return false;
        }
    }
    read_num_t read_id;
    read_length_t start;
};

struct Terminal_node_info {
    Terminal_node_info(): node_type(NO_PAIR), align_read_id(IMPOSSIBLE_READ_NUM){}
    Terminal_node_info(uint8_t node_type_, read_num_t align_read_id_) :
                           node_type(node_type_), align_read_id(align_read_id_) {}
    uint8_t node_type;
    read_num_t align_read_id; // two read

};

//search info
struct S_info {
    S_info() : color(WHITE), parent(NULL), d(0), mid_seq_len(0) {}
    uint8_t color;
    void * parent;
    int d;
    int mid_seq_len;
};

struct bundled_node_p {
    bundled_node_p(const std::string & seq_="",
                   const node_id_t & node_id_=IMPOSSIBLE_NODE_ID):
                   seq(seq_), node_id(node_id_), count(1), norm{1},
                   copy_count(0) {}

    inline read_length_t seq_len(){return seq.size();}
    inline std::string to_string()
    {
        return seq + "\t" + std::to_string(count);
    }
    node_id_t node_id;
    std::string seq;
    kmer_count_node_t count;
    kmer_count_node_t prevalence;
    kmer_count_node_t copy_count;
    kmer_count_node_t norm;
    
    std::vector<Bdg_read_info> reads_info;
    std::vector<Terminal_node_info> term_nodes_info;
    S_info s_info;
};

struct bundled_edge_p {
    bundled_edge_p(const edge_weight_t & weight_=IMPOSSIBLE_EDGE_WEIGHT_NUM,
                   const double & count_=0)
                        : weight(weight_), count(count_) {}
    edge_weight_t weight; // is overlap
    double count;   // is the count
};

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


bool acyclic_check();
bool DFS_visit(vd_t);

#endif	/* SEQ_GRAGH_TYPE_H */
