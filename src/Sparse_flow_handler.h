/*
 * File:   Sparse_flow_handler.h
 * Author: bx
 *
 * Created on November 3, 2017, 12:14 PM
 */

#ifndef SPARSE_FLOW_HANDLER_H
#define	SPARSE_FLOW_HANDLER_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/any.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <functional>
#include <numeric>
#include <string.h>
#include <fstream>
#include <glpk.h>
#include <random>
#include <stdlib.h>
#include <pthread.h>

#include "shc_type.h"
#include "json_parser.h"
#include "local_file_structure.h"
#include "shc_google_sparsehash.h"
//#include "Sequence_graph_handler.h"
#include "log.h"
#include "Seq_gragh_type.h"

#define SF_WORK_STEP_UPDATE 30

//debug
//#define LOG_LP_SUMMARY
//#define LOG_FLOW
//#define LOG_SF

#define EMPTY_NODE_ID (std::pow(2,sizeof(node_id_t)*8-1)-1)
#define DELETE_NODE_ID (std::pow(2,sizeof(node_id_t)*8-2)-1)
#define IGNORE 0.0

struct Comp_graph {
    Comp_graph(int comp_i_, int graph_i_): comp_i(comp_i_), graph_i(graph_i_) {}
    int comp_i;
    int graph_i;
};

class Sparse_flow_handler;
struct SF_work_info {
    SF_work_info(std::deque<Comp_graph> * work_list_,
                 Sparse_flow_handler * context_):
                 context(context_), work_list(work_list_) {}
    std::deque<Comp_graph> * work_list;
    Sparse_flow_handler * context;
    pthread_mutex_t * work_lock_ptr;
    pthread_mutex_t * write_lock_ptr;
    int process_i;
};

struct Sparse_flow_setting {
    int path_sparsity;
    bool use_Y_path;
    int num_trials;
    double eps;
    double removal_factor;
    double sparsity_factor;
    double tol;
    bool is_use_topo_reverse;
    int multiple_test;

    std::string output_path;
};

double norm_1_helper(double x, double y);
template<class T1, class T2> double dot_product(T1 & v1, T2 & v2);

class Sparse_flow_handler {
public:
    struct Path_info {
        Path_info(int id_, int location_) : id(id_), location(location_) {}
        int id;
        int location; //start from 0, unit is number of vd
    };

    typedef std::pair<int, double> index_flow_t;
    struct Index_flow_sorter {
        bool operator()(const index_flow_t& p1, const index_flow_t& p2) const {
    	return( p1.second > p2.second);
        }
    };


    typedef google::dense_hash_map<node_id_t, vd_t, hash_u64, equ64> Node_id_vd_map;
    typedef google::dense_hash_map<node_id_t, vd_t, hash_u64, equ64>::iterator Node_id_vd_map_iterator;

    typedef std::multimap<vd_t, Path_info> Node_path_MMap;
    typedef std::multimap<vd_t, Path_info>::iterator Node_path_MMap_iterator;
    typedef std::pair<Node_path_MMap_iterator,Node_path_MMap_iterator>
                                                Node_path_MMap_iterator_pair;

    typedef std::map<vd_t, vd_t> Vd_Vd_map;
    typedef std::map<vd_t, vd_t>::iterator Vd_Vd_map_iterator;

    typedef std::map<vd_t, std::deque<vd_t> > Vd_Vds_map;
    typedef std::map<vd_t, std::deque<vd_t> >::iterator Vd_Vds_map_iterator;

    Sparse_flow_handler(Shannon_C_setting & setting_):
                        setting(setting_), distribution(0,1.0)
    {
        nodeId_vd_map.set_deleted_key(DELETE_NODE_ID);
        nodeId_vd_map.set_empty_key(EMPTY_NODE_ID);
        sf_setting.use_Y_path = true;

        glob_node_id = 0;


        sf_setting.eps = 0.001;
        sf_setting.path_sparsity = 10;
        sf_setting.removal_factor = 0.4;
        sf_setting.sparsity_factor = 0.4;
        //sf_setting.is_use_topo_reverse = is_use_reverse_;
        sf_setting.multiple_test = setting.sparse_flow_setup.multiple_test;
        //generator.seed(31);
        //modify_writer.open(setting.local_files.output_path+"/modify.log");
    }

    //used for multi-thread
    void * run_sparse_flow_handler_helper( Comp_graph comp_graph,
            std::string node_path, std::string edge_path, std::string path_path,
            pthread_mutex_t * write_lock_ptr, std::string output_path);

    // to be implemented
    void process_all_graph(int comp_i);
    void process_one_graph_component(int graph_i, std::string & node_path,
                              std::string & edge_path, std::string & path_path);
    void load_node_edge_path_files();
    void clear();
    bool acyclic_check();
    int get_num_nodes() {return boost::num_vertices(graph);}
    int get_num_edges() {return boost::num_edges(graph);}

    void log_all_node(bool is_log_node);
    void log_node_info(vd_t vd, bool is_log_nodes);

    void log_flow(std::vector<double> & flow, std::vector<vd_t> & in_nodes,
                                              std::vector<vd_t> & out_nodes);


private:
    // reconstructing graph
    bool load_graph(std::string & node_file, std::string & edge_file);
    bool load_node(std::string & node_file);
    bool load_edge(std::string & edge_file);
    void load_path(std::string & path_file);

    bool load_py_graph(std::string & node_file, std::string & edge_file);
    bool load_py_node(std::string & node_file);
    bool load_py_edge(std::string & edge_file);
    void load_py_path(std::string & path_file);

    // add pseudo source, sink nodes
    void add_start_end_node();

    void link_all_Y_path(std::set<vd_t> & remove_vd);
    bool link_Y_path(vd_t vd);
    bool is_connect_start_node(vd_t vd);
    bool is_connect_end_node(vd_t vd);
    // sparse flow
    void sparse_flow_on_all_nodes();
    bool sparse_flow_on_one_node(vd_t vd, std::vector<double> & flow,
            std::vector<double> & in_counts, std::vector<double> & out_counts,
            std::vector<int> & sub_flow_indicator, std::vector<vd_t> & in_nodes,
                                                      std::vector<vd_t> & out_nodes);
    // get known path indicator
    vd_t get_ori_vd(vd_t vd_source);
    void get_sub_flow_indicator(vd_t vd, std::vector<int> & sub_flow_indicator,
                   std::vector<vd_t> & in_nodes, std::vector<vd_t> & out_nodes);
    void make_in_out_flow_equal(std::vector<double> & in_counts,
                                             std::vector<double> & out_counts);
    double path_decompose(vd_t vd, std::vector<double>& gamma,
                          std::vector<double> & equ_constraint,
            std::vector<int>& sub_flow_indicator, std::vector<double> & flow, double & scale);
    void truncate_flow_value(int num_in, int num_out,
             std::vector<double> & flow, std::vector<double> & flow_2,
             std::vector<double> & in_counts, std::vector<double> & out_counts);

    void old_truncate_flow_value(int num_in, int num_out,
             std::vector<double> & flow, std::vector<double> & flow_2,
             std::vector<double> & in_counts, std::vector<double> & out_counts);

    void py_make_flow_sparser(int product, std::vector<double> & flow);

    bool is_update_flow(std::vector<double> & flow, std::vector<double> & min_flow,
              std::vector<int>& sub_flow_indicator, double min_unknown_flow);
    // modify graph after sparse flow
    void modify_local_graph(vd_t vd, std::vector<double> & flow,
               std::vector<vd_t> & in_nodes, std::vector<vd_t> & out_nodes);
    void link_leaf_node_to_terminal(std::vector<vd_t> & in_nodes,
                                                std::vector<vd_t> & out_nodes);
    void link_node_to_terminals(vd_t vd);
    void condense_nodes(std::vector<vd_t> & nodes, std::set<vd_t> & remove_vd);

    void update_sorted_vds(bool is_reverse);
    void remove_decomposed_node(std::set<vd_t> & remove_vd);
    // output path
    void dump_SF_seq(std::string & out_file);
    void output_seq_helper(std::ofstream & file_writer, std::string & curr_seq,
                    std::vector<int> & delimit, std::vector<vd_t> & stack_vd);

    //multi pthread
    void dump_SF_seq(std::string & out_file, pthread_mutex_t * writer_lock_ptr);
    void output_seq_helper_mt(std::ofstream & file_writer, pthread_mutex_t * writer_lock_ptr,
                     std::string & curr_seq, std::vector<int> & delimit,
                     std::vector<vd_t> & stack_vd);
    // check acyclic
    bool DFS_visit(vd_t);

    void update_contain_map(vd_t vd, vd_t vd_contained);
    bool is_subflow_on_vd_unknown(vd_t vd, vd_t u, vd_t w, int b_d, int a_d);


    //elementary helper function
    void get_in_out_edges(vd_t vd,
             std::vector<vd_t> & in_nodes, std::vector<vd_t> & out_nodes,
             std::vector<double> & in_counts, std::vector<double> & out_counts);
    double get_norm_2(std::vector<double> & v);
    double get_norm_1(std::vector<double> & v);
    int get_norm_0(std::vector<double> & v);
    int get_num_nonzero_flow(std::vector<double> & flow,
                                    std::vector<int> & sub_flow_indicator);
    int get_num_nonzero_flow(std::vector<double> & flow);
    double get_flow_difference(std::vector<double> & flow1,
                               std::vector<double> & flow2);
    void get_i_j(int flow_index, int num_out, int & i, int & j);
    void set_trials(int m, int n);
    double sum_vector(std::vector<double> & v);
    int sum_vector(std::vector<int> & v);
    double get_scale(std::vector<double> & in_counts, std::vector<double> & out_counts);

    //debug
    void log_graph_struct (bool show_seq);
    void log_term_node(bool is_front, bool is_back, bool detail);
    void log_local_node(vd_t vd, bool is_show_seq);
    void log_around_local_node(vd_t vd,
                    std::vector<vd_t> & in_nodes, std::vector<vd_t> & out_nodes);
    void check_any_condensible_node();

    vd_t local_condense(vd_t vd, std::set<vd_t> & vd_remove_set);
    void condense_left(vd_t vd, vd_t vd_new, std::set<vd_t> & vd_remove_set);
    void condense_right(vd_t vd, vd_t vd_new, std::set<vd_t> & vd_remove_set);
    bool condense(vd_t vd, bool us_check_left , vd_t vd_new);
    bool is_condensable(vd_t vd, bool is_left);
    void condense_graph(std::set<vd_t> & vd_remove_set);


    seq_graph_t graph;
    std::vector<vd_t> sorted_vds;
    Node_id_vd_map nodeId_vd_map;
    std::vector< std::vector<vd_t> > known_paths;
    Node_path_MMap paths_for_node;

    //Vd_Vd_map new_to_ori;
    Vd_Vds_map vd_contains_map;

    vd_t vd_start;
    vd_t vd_end;
    node_id_t glob_node_id;

    int comp_id;
    int graph_id;
    int num_out_seq;

    Sparse_flow_setting sf_setting;
    Shannon_C_setting & setting;

    std::normal_distribution<double> distribution;
    std::default_random_engine generator;

    //std::ofstream modify_writer;
};

struct Sparse_flow_works {
    Sparse_flow_works () {}
    Sparse_flow_works (std::deque<Comp_graph> * work_list_, Shannon_C_setting & setting_,
                    pthread_mutex_t * work_lock_ptr_) :
            work_list(work_list_),  setting(setting_),
            work_lock_ptr(work_lock_ptr_) {}
    std::deque<Comp_graph> * work_list;
    Shannon_C_setting setting;
    pthread_mutex_t * work_lock_ptr;
    pthread_mutex_t * write_lock_ptr;
    int process_i;
    std::string output_path;
    int init_total_work;

    static void * run_multi_sparse_flow(void * sparse_flow_work_ptr)
    {
        Sparse_flow_works sparse_flow_work =
            (*static_cast<Sparse_flow_works * >(sparse_flow_work_ptr));
        Sparse_flow_handler * sparse_flow_ptr =
                new Sparse_flow_handler(sparse_flow_work.setting);
        while(1)
        {
            // take the work
            pthread_mutex_lock(sparse_flow_work.work_lock_ptr);
            if(sparse_flow_work.work_list->empty())
            {
                //std::cout <<"releasing a thread " << std::endl;
                pthread_mutex_unlock(sparse_flow_work.work_lock_ptr);
                return ((void*) 0);
            }
            Comp_graph comp_graph = sparse_flow_work.work_list->front();
            sparse_flow_work.work_list->pop_front();
            int num_comp_left = sparse_flow_work.work_list->size();
            pthread_mutex_unlock(sparse_flow_work.work_lock_ptr);

            if(num_comp_left%SF_WORK_STEP_UPDATE==0)
            {
                int percentage = (100 * (sparse_flow_work.init_total_work-num_comp_left))/sparse_flow_work.init_total_work;
                printf("[ %3d%%] SF process %4d, \t num left %d/%d\n",
                        percentage, sparse_flow_work.process_i, num_comp_left,
                        sparse_flow_work.init_total_work);
            }

            int comp_id = comp_graph.comp_i;
            int graph_id = comp_graph.graph_i;
            // prepare files for process
            Local_files & lf = sparse_flow_work.setting.local_files;
            std::string node_path(lf.output_seq_graph_path +
                                         lf.node_prefix + std::to_string(comp_id) +
                                         "/node" + std::to_string(graph_id));
            std::string edge_path(lf.output_seq_graph_path +
                                         lf.edge_prefix + std::to_string(comp_id) +
                                         "/edge" + std::to_string(graph_id));
            std::string path_path(lf.output_seq_graph_path +
                                         lf.path_prefix + std::to_string(comp_id) +
                                         "/path" + std::to_string(graph_id));

            sparse_flow_ptr->run_sparse_flow_handler_helper( comp_graph,
                node_path, edge_path, path_path,
                sparse_flow_work.write_lock_ptr, sparse_flow_work.output_path);
        }

        return NULL;
    }
};


#endif	/* SPARSE_FLOW_HANDLER_H */
