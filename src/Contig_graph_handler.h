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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include "local_file_structure.h"
#include <memory>
#include <stdlib.h>


#include "shc_type.h"

#include "encoding.h"

#include "Contig_handler.h"
#include "Kmer_handler.h"
#include "shc_google_sparsehash.h"
#include <metis.h>
#include "json_parser.h"

//#include "Fasta_entry.h"
//#include "Fasta_reader.h"

#define ONLY_EDGE_WEIGHT "001"
#define METIS_NODE_OFFSET 1
#define EMPTY_KEY (std::pow(2,sizeof(contig_num_t)*8-1)-1)
#define METIS_IMBALANCE_PRECISION 1000.0

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

    typedef google::dense_hash_map<contig_num_t, vd_t, hash_contig, eq_contig> Contig_vertex_map;
    typedef google::dense_hash_map<contig_num_t, vd_t, hash_contig, eq_contig>::iterator Contig_vertex_map_iterator;
    typedef google::dense_hash_map<contig_num_t, vd_t, hash_contig, eq_contig>::const_iterator Contig_vertex_map_const_iterator;

    struct Metis_input {
        Metis_input() {xadj.push_back(0);}
        void metis_reset_data()
        {
            adjncy.clear();
            xadj.clear();
            adjwgt.clear();
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

public:
    Contig_graph_handler(kmer_len_t length, Kmer_handler * khp,
         Contig_handler * chp, Local_files *lfp, Metis_setup & ms):
         k1mer_len(length), kh(khp), ch(chp),
         component_array(chp->num_contig, IMPOSSIBLE_COMP_NUM), metis_setup(ms)
         {
             contig_vertex_map.set_empty_key(EMPTY_KEY);
             explorable_contig_set.set_empty_key(ch->num_contig+2);
             explorable_contig_set.set_deleted_key(ch->num_contig+3);
             curr_component_num = 0;
             is_set_collect_comp_num = false;
             lf = lfp;
             if(metis_setup.is_multiple_partition)
             {
                 component_array_aux.resize(chp->num_contig);
                 memset(&component_array_aux.at(0), IMPOSSIBLE_COMP_NUM, chp->num_contig);
             }
         };

    void remove_read_duplicate(std::string out_dir_path);

    void group_components();
    void break_and_keep_component();
    void run_single_pass_metis(std::vector<comp_num_t> & array, std::vector<idx_t> * part);
    void add_weight_to_cut_and_update_graph(std::vector<idx_t> * part);


    void get_connected_contig_index(Kmer_counter_map_iterator & it);
    void increment_edge_weight(contig_num_t i, contig_num_t j);

    void create_metis_array_input();
    void dump_component_array(std::string & filename);

    void dump_graph_to_metis_file();
    void create_metis_file_format_from_graph(std::string & filename);

    void assign_reads_to_components(int num_test, std::string input_read_path);
    void assign_kmer_to_components();

    void assign_paired_read_to_components(int num_test);

    //void merge_group();
    //debug
    void log_metis_input_data();


private:
    graph_t graph;
    kmer_len_t k1mer_len;
    Contig_vertex_map contig_vertex_map;
    std::vector<char> curr_contig;

    std::stack<contig_num_t> contig_stack;
    Contig_set explorable_contig_set;

    std::vector<comp_num_t> component_array;
    std::vector<comp_num_t> component_array_aux;

    comp_num_t curr_component_num; // after processing the largest comp num is this number -1

    comp_num_t collect_comp_num;
    bool is_set_collect_comp_num;

    Kmer_handler * kh;
    Contig_handler * ch;

    struct Metis_input metis_input;
    struct Metis_setup metis_setup;

    Local_files * lf;
};

#endif	/* CONTIG_GRAPH_HANDLER_H */
