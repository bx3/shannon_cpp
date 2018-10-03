/*
 * File:   Sequence_graph_handler.h
 * Author: bx
 *
 * Created on October 2, 2017, 8:15 PM
 */
#ifndef SEQUENCE_GRAPH_HANDLER_H
#define	SEQUENCE_GRAPH_HANDLER_H

#include <pthread.h>
#include <boost/algorithm/string/trim.hpp>
#include <vector>
#include <string.h>
#include <iostream>
#include <fstream>
#include <set>
#include <deque>
#include <stack>
#include <queue>
#include <list>

#include "shc_google_sparsehash.h"
#include "shc_type.h"
#include "json_parser.h"
#include "Kmer_handler.h"
#include "local_file_structure.h"
#include "log.h"
#include "Collect_reads.h"
#include <cmath>
#include "Seq_gragh_type.h"
#include "Read_subsampler.h"
//#include "Multi_graph_handler.h"

#define EMPTY_READ_NUM (std::pow(2, sizeof(read_num_t)*8)-1)
#define EMPTY_READ_LENGTH (std::pow(2, sizeof(read_length_t)*8)-1)
#define TWO_POWER_16 65536
#define TWO_POWER_8 256
#define POINTER_HASH 2654435761
#define IMPOSSIBLE_READ_HASH (std::pow(2, sizeof(uint64_t)*8-2)-1)

#define PAIR_1 1
#define PAIR_2 2

#define KNOWN_PATH_VEC_INIT_SIZE 10

#define SEQ_GRAPH_WORK_STEP_UPDATE 5
#define DEFAULT_COMP_MEM -2
#define DEFAULT_COMP_NUM -2
#define MEM_NOT_ENOUGH_ID -3

//#define SHOW_SEQ_GRAPH_PROGRESS
//#define LOG_SEQ_GRAPH
//#define LOG_PAIR_SEARCH

//#define PRINT_GRAPH_STATS
//#define PRINT_TIME

class Sequence_graph_handler;

struct Work_info {
    Work_info () {}
    ~Work_info () {}
    Work_info (std::deque<int> * work_list_, Sequence_graph_handler * context_):
                        work_list(work_list_), context(context_) {}
    std::deque<int> * work_list;
    Sequence_graph_handler * context;
    pthread_mutex_t * work_lock_ptr;
};

struct Special_seq_input {
    std::string kmer_path;
    std::string read_path_single_prefix;
    std::string read_path_p1_prefix;
    std::string read_path_p2_prefix;
};


class Sequence_graph_handler {
    friend class Kmer_handler;
    friend class Multi_graph_handler;

public:
    //Google sparsehash
    typedef google::dense_hash_map<std::string, vd_t,
                                std::hash<std::string>, eqstr> Kmer_Node_map;
    typedef google::dense_hash_map<std::string, vd_t,
             std::hash<std::string>, eqstr>::iterator Kmer_Node_map_iterator;

    typedef google::dense_hash_map<read_num_t, std::vector<vd_t>,
                                    hash_u32, equ32> Read_path_map;
    typedef google::dense_hash_map<read_num_t, std::vector<vd_t>,
                        hash_u32, equ32>::iterator Read_path_map_iterator;

    typedef tsl::hopscotch_set<vd_t> Vd_set;
    typedef tsl::hopscotch_set<vd_t>::iterator Vd_set_iterator;
    typedef tsl::hopscotch_set<vd_t>::const_iterator Vd_set_const_iterator;

    struct hash_read_info
    {
        uint64_t operator() (const Bdg_read_info& read_info) const
        {
            return(static_cast<uint64_t>(read_info.read_id)<<31 | static_cast<uint64_t>(read_info.start));
        }
    };

    struct equ_read_info
    {
        bool operator()(const Bdg_read_info & read_info1, const Bdg_read_info & read_info2) const
        {
            return (read_info1.read_id==read_info2.read_id && read_info1.start==read_info2.start);
        }
    };

    typedef google::dense_hash_set<Bdg_read_info, hash_read_info, equ_read_info> Read_set;
    typedef google::dense_hash_set<Bdg_read_info, hash_read_info, equ_read_info>::iterator Read_set_iterator;

    typedef google::dense_hash_set<read_num_t, hash_read, eq_read> Read_unique_set;
    typedef google::dense_hash_set<read_num_t, hash_read, eq_read>::iterator Read_unique_set_iterator;

    typedef std::set<vd_t> Node_set_t;
    typedef std::set<vd_t>::iterator Node_set_iterator;


    typedef std::pair<vd_t, kmer_count_node_t> Node_Count_Pair;
    struct Node_count_sorter {
        bool operator()(const Node_Count_Pair& p1, const Node_Count_Pair& p2) const {
            return( p1.second > p2.second);
        }
    };

    struct Node_index {
        Node_index(vd_t vd_, read_length_t start_) : vd(vd_), start(start_) {}
        vd_t vd;
        read_length_t start;
    };

    struct Node_pair {
        Node_pair(vd_t vd_start_, vd_t vd_stop_):
                       vd_start(vd_start_), vd_stop(vd_stop_){}
        bool operator < (const Node_pair& np) const
        {
            return vd_start < np.vd_start ||
                   (!(np.vd_start<vd_start) && (vd_stop < np.vd_stop));
        }
        vd_t vd_start;
        vd_t vd_stop;
    };
    struct equ_node_pair {
        bool operator()(const Node_pair& np1, const Node_pair& np2) const {
            return( np1.vd_start==np2.vd_start && np1.vd_stop==np2.vd_stop);
        }
    };
    struct hash_node_pair {
        size_t operator() (const Node_pair& np) const {
            return((size_t)np.vd_start); //*POINTER_HASH + (size_t)np.vd_stop
        }
    };

    struct Node_struct_info {
        Node_struct_info(int in_d_, int out_d_) : in_d(in_d_), out_d(out_d_) {}
        bool operator < (const Node_struct_info& nsi) const
        {
            return std::make_pair(in_d, out_d) < std::make_pair(nsi.in_d, nsi.out_d);
        }
        int in_d;
        int out_d;
    };

    typedef std::map<Node_pair, std::vector<vd_t>,
             equ_node_pair> Node_pair_path_map;
    typedef std::map<Node_pair, std::vector<vd_t>,
             equ_node_pair>::iterator Node_pair_path_map_iterator;

    Sequence_graph_handler(Shannon_C_setting & setting_) :
                setting(setting_),
                coll_read_list(setting_.has_pair, setting_.has_single,
                setting_.single_read_length, setting_.pair_1_read_length,
                setting_.pair_2_read_length, setting_.is_compress, setting_)
    {
        has_single = setting_.has_single;
        has_pair = setting_.has_pair;
        kmer_length = setting_.kmer_length;
        max_hop_pair = setting_.seq_graph_setup.max_hop_pair;
        max_hop_path = setting_.seq_graph_setup.max_hop_path;
        special_input = false;


        glob_node_id = NODE_ID_NORMAL_START;
        is_multi_thread = false;
        if(has_single && !has_pair)
            size_threshold = setting_.single_read_length;
        else if(!has_single && has_pair)
            size_threshold =
                  (setting_.pair_1_read_length + setting_.pair_2_read_length)/2;
        else if(has_single && has_pair)
            size_threshold = (setting_.pair_1_read_length +
                    setting_.pair_2_read_length + setting_.single_read_length)/3;
        else
            exit(1);


        prevalence_threshold = 1;
        hamming_frac = 0.1;
        peak_mem = 0;

        read_empty = Bdg_read_info(IMPOSSIBLE_READ_NUM, IMPOSSIBLE_READ_LENGTH);
        read_delete = Bdg_read_info(IMPOSSIBLE_READ_NUM-1000, IMPOSSIBLE_READ_LENGTH-1000);

        //debug_file
        num_read_on_node = 0;
        num_condensed_rm = 0;
        total_num_bridged =0;

        min_read_len = IMPOSSIBLE_READ_LENGTH;
        if(has_single)
        {
            min_read_len = setting.single_read_length;
        }
        if(has_pair)
        {
            min_read_len = std::min(min_read_len, setting.pair_1_read_length);
            min_read_len = std::min(min_read_len, setting.pair_2_read_length);
        }

        //bridge_nodes_writer.open(setting.local_files.output_path+"/bridge_node.log");
        //graph_writer.open(setting.local_files.output_path+"/graph.log");
        //reads_writer.open(setting.local_files.output_path+"/reads.log");
        //cycle_writer.open(setting.local_files.output_path+"/cycle.log");

        //rm_sus_writer.open(setting.local_files.output_path+"/sus.log");
    }

    // available public function
    int64_t run_it(int comp_i, bool is_single_component);
    void setup_input_file_path(int comp_i);
    void setup_input_file_path(std::string kmer_path, std::string s_read_path,
                            std::string p1_read_path, std::string p2_read_path);
    void build_kmer_graph_from_reads();
    void build_kmer_graph_from_edges();
    void load_all_read(Kmer_Node_map & kmer_node_map);
    void condense_graph();
    void assign_read_to_xnode();
    void update_xnode_set();
    void bridge_all_xnodes();
    void find_approximate_copy_count();
    void find_copy_count();
    void break_self_loops();
    void break_all_cycles();
    bool acyclic_check();
    int find_known_path(int curr_hop);
    uint64_t resolve_all_pair_reads(int search_hop, std::set<vd_t> & check_vd);
    void output_components(std::string & node_output_path,
                           std::string & edge_output_path,
                           std::string & path_output_path);
    void collapse_all();
    void remove_all_suspicious_nodes();
    uint64_t find_all_simple_pair();
    void get_all_descendant(std::set<vd_t> & repeat_set, vd_t root_vd);

    std::deque<vd_t>::iterator thread_safe_find(
                         std::deque<vd_t> & cycle_path, vd_t vd_target);


    int get_num_nodes() {return boost::num_vertices(graph);}
    int get_num_edges() {return boost::num_edges(graph);}
    int get_num_isolated_nodes();
    int get_num_xnodes();
    void set_curr_comp(int curr_comp_)
    {
        curr_comp = curr_comp_;
        coll_read_list.set_comp_id(curr_comp_);
    }

    // for debug
    void log_all_node(bool is_log_node, bool is_log_read);
    void log_all_edge(bool is_log_weight, bool is_log_count);
    void log_xnode_set(bool is_log_node, bool is_log_read);
    void log_read_seq(Read_acc acc);
    void log_node_info(vd_t vd, bool is_log_nodes , bool is_log_read_info, bool is_log_bfs);
    void log_seq_only(vd_t vd, read_length_t start);
    void log_node_bfs_info(vd_t vd, bool is_show_color, bool is_shown_depth);
    void log_edge_count();
    void log_term_array(bool is_show_seq);
    void log_graph_to_file(std::ofstream & writer, int id);
    void log_all_nodes_reads_to_file(std::ofstream & writer);

    void log_classify_edge_types();
    void log_classify_node_types(std::ofstream & writer, int id);


private:

    //load reads
    void load_all_single_read(std::string& read_path, Kmer_Node_map & kmer_node_map);
    void load_all_paired_read(Kmer_Node_map & kmer_node_map);
    void load_all_paired_read_no_concat(std::string read_path_p1,
                    std::string read_path_p2, Kmer_Node_map & kmer_node_map);
    vd_t update_edge_count_with_read(std::string & base, Kmer_Node_map & kmer_node_map);
    vd_t traverse_read_is_all_kmer_node_valid(std::string & base,
                              Kmer_Node_map & kmer_node_map, uint8_t node_type);

    // specify which graph to work on
    void clear_seq_graph_mem();
    // build kmer graph
    void use_reads_to_build_edge(std::string& read_path,
                    Kmer_Node_map& kmer_node_map, read_num_t *read_id);
    void use_pair_reads_to_build_edge(std::string& read_path_1, std::string& read_path_2,
         Kmer_Node_map& kmer_node_map, read_num_t *read_id);

    vd_t process_one_read_to_build_graph( Kmer_Node_map& kmer_node_map,
                                std::string & read_base, uint8_t node_type);

    //remove suspicous
    bool remove_suspicious_nodes(uint64_t & node_removed);
    bool is_suspicious(vd_t vd);
   // collapse
    bool collapse_nodes(uint64_t & node_removed);
    bool collapse_successor(vd_t vd, std::set<vd_t> & vd_remove);
    vd_t collapse(vd_t vd_1, vd_t vd_2);
    bool similar(vd_t vd_1, vd_t vd_2);

    // for searching unique path between pair read
    bool is_node_contain_terminals(vd_t vd);
    void clear_term_info(uint64_t align_index);
    bool check_and_add_simple_read_pair(uint64_t i, bool & is_sr);
    bool find_middle_seq(vd_t vd, read_num_t local_read_1, read_num_t local_read_2,
                std::string & middle, int64_t & middle_len, bool & is_mismatch);
    int64_t find_terminal_index(vd_t vd, std::string & vd_seq,
            char * read_ptr, read_length_t len, read_length_t acc_len, bool is_p1);

    bool search_path_in_pair_nodes_with_bfs(vd_t vd1, vd_t vd2, int max_hop,
                              std::vector<vd_t> & path, Node_pair_path_map & node_pair_path_map, read_num_t l_start_read_id);
    bool is_unique_path_to_u(vd_t u, vd_t & repeat_vd);

    void paint_queue_node(vd_t vd_root, vd_t vd_check, std::deque<vd_t> & vd_queue);
    bool is_ancestrial(vd_t vd_root, vd_t vd_stop, vd_t vd_start);
    void create_path_vd_set(vd_t vd_end, vd_t vd_start, std::set<vd_t> & vd_set);
    bool check_and_create_path_vector(vd_t vd_end, vd_t vd_start, std::vector<vd_t> & vd_list);
    void convert_vd_path_to_string(std::vector<vd_t> & path, std::string & str,
                              read_num_t local_read_1, read_num_t local_read_2);
    bool search_next_layer(vd_t vd1, vd_t vd2, bool is_forward);
    // bridge xnodes
    void get_xnodes(std::set<vd_t> & new_nodes, std::set<vd_t> & xnodes);
    bool is_xnode_bridged(vd_t vd);
    void perform_bridging(vd_t vd);
    bool is_read_bridge_node(read_num_t read_id, read_length_t start, vd_t vd);
    vd_t create_left_node(ed_t ed);
    vd_t create_right_node(ed_t ed);
    // condense simple
    vd_t merge_two_vertices(vd_t vd1, vd_t vd_2);
    void transfer_read_align_index(vd_t vd1, vd_t vd);
    void erase_terminal_node_info(vd_t vd);
    // local condense
    void local_condense_bridged(std::vector<vd_t> & u_list, std::vector<vd_t> & w_list);
    vd_t local_condense(vd_t vd, std::set<vd_t> & vd_remove_set);
    void condense_left(vd_t vd, vd_t vd_new, std::set<vd_t> & vd_remove_set, Read_set & read_set);
    void condense_right(vd_t vd, vd_t vd_new, std::set<vd_t> & vd_remove_set, Read_set & read_set);
    bool condense(vd_t vd, bool us_check_left , vd_t vd_new,
              Read_set & read_set, int & seq_len_before_last_concat);
    bool is_condensable(vd_t vd, bool is_left);
    bool py_special_condense(vd_t vd);

    // path search
    bool is_partial_match(const char * ptr1, const read_length_t len1,
                           const char * ptr2, const read_length_t len2);
    bool search_seq(const char * seq_ptr, read_length_t seq_len,
                     vd_t vd, read_length_t start,
                     int curr_hop, std::vector<vd_t> & path);
    //break cycles
    bool find_cycle(std::set<vd_t> & acyclic_node_set, std::deque<vd_t> & cycle_path);
    bool is_node_inside_cycle(vd_t vd, std::set<vd_t> & acyclic_node_set, std::deque<vd_t> & cycle_path);
    void simple_break_cycle(std::deque<vd_t> & cycle_path);
    bool break_cycle(std::deque<vd_t> cycle_path);
    void break_cycle_splitting(vd_t vd_split, vd_t vd_c_in, vd_t vd_c_out,
        int sum_copy_count, int sum_c_in_and_o_out);
    bool DFS_visit(vd_t);

    void py_break_cycles();
    bool py_find_cycle(std::set<vd_t> & no_cycle, std::deque<vd_t> & traversed);
    void py_reachable_cycle(vd_t vd, std::set<vd_t> & no_cycle,
                std::vector<vd_t> & traversed, std::deque<vd_t> & cycle);

    // output dump
    void seq_graph_output_dir_setup_helper(std::string & dir);
    void add_component(vd_t vd, std::set<vd_t> & nodes, std::set<ed_t> & edges);
    bool variant_topological_sort(std::set<vd_t> & nodes, std::vector<vd_t> & sorted_nodes);
    bool is_all_predecessors_included(vd_t vd_root, std::set<vd_t> & added_nodes);
    void path_to_string(std::vector<vd_t> & path, std::string & out_string);
    std::string edge_to_string(ed_t ed);
    // preprocess
    void pre_process_read();


    void remove_nodes_edges_if_not_cover_by_reads();

    // elementary help function
    void refresh_bridging_reads(vd_t vd);
    bool is_read_bridge_left(Bdg_read_info & read_info, vd_t vd);
    bool is_read_bridge_right(Bdg_read_info & read_info, vd_t vd);
    void link_in_edges(vd_t vd1, vd_t vd);
    void link_out_edges(vd_t vd2, vd_t vd);
    bool is_xnode(vd_t vd);
    bool is_read_bridge_node(std::string & read_seq, read_length_t start, vd_t vd);
    bool is_read_bridge_node(const char * read_ptr, read_length_t len,
                                        read_length_t info_start, vd_t vd);
    inline void node_link_read(vd_t vd, read_num_t read_id, read_length_t start)
    {
        graph[vd].reads_info.emplace_back(read_id, start);
    }
    void get_successors(vd_t vd, std::set<vd_t> & successors);
    void get_predecessors(vd_t vd, std::set<vd_t> & predecessors);
    void remove_set_nodes(std::set<vd_t> & vd_remove);
    bool trim_read(std::string & read, uint64_t i_read);
    void reset_node_color(std::set<vd_t> & colored_set);
    bool is_only_self_loop(vd_t vd);
    void special_condense_self_loop();
    bool is_contain_self_loop(vd_t vd, ed_t & self_edge);
    size_t get_total_num_reads(std::string file_path_prefix);

    inline double average_prevalence(vd_t vd)
    {
        return static_cast<double>(graph[vd].prevalence) / graph[vd].count;
    }
    //not used
    void load_kmer_and_build_graph();


    // should be cleared for each component
    seq_graph_t graph;
    struct Collect_reads coll_read_list;

    int curr_comp;

    std::set<vd_t> xnode_set;
    std::set<vd_t> paired_node_set;

    //specify the node which serves as i th read terminal reads
    std::vector<vd_t> p1_end;
    std::vector<vd_t> p2_front;

    //store known path
    std::set<std::vector<vd_t> > known_path_set;
    std::map<ed_t, uint64_t> known_edges;

    // should be updated
    node_id_t glob_node_id;

    //settings
    std::string kmer_path;
    std::string read_path_single_prefix;
    std::string read_path_p1_prefix;
    std::string read_path_p2_prefix;

    Block_timer timer;

    uint8_t kmer_length;
    Shannon_C_setting & setting;
    int max_hop_path;
    int max_hop_pair;
    bool has_pair;
    bool has_single;

    int size_threshold;
    int prevalence_threshold;
    float hamming_frac;

    bool is_multi_thread;
    bool special_input;

    int working_component;

    Bdg_read_info read_empty;
    Bdg_read_info read_delete;

    bool is_run_single_component;

    //debug
    uint64_t curr_align_index;
    std::set<vd_t> exist_vd;
    Block_timer test_timer;

    uint64_t num_read_on_node;
    uint64_t num_condensed_rm;
    Block_timer mb_condense_timer;
    uint64_t total_num_bridged;
    uint64_t num_single_seq;

    std::ofstream bridge_nodes_writer;
    std::ofstream graph_writer;
    std::ofstream reads_writer;
    std::ofstream cycle_writer;
    std::ofstream nodes_struct_writer;

    //std::ofstream rm_sus_writer;
    Block_timer resolve_pair_timer;

    uint64_t num_edge_added;
    uint64_t num_assigned_read;

    read_length_t min_read_len;

    int64_t peak_mem;
};


struct Seq_graph_works {
    Seq_graph_works ()
    {
        finish_token = "-1";
    }
    Seq_graph_works (std::deque<int> * work_list_, Shannon_C_setting & setting_,
                    pthread_mutex_t * work_lock_ptr_,
                    int read_fd_, int wk_fd_) :
            work_list(work_list_),  setting(setting_),
            work_lock_ptr(work_lock_ptr_),
            read_fd(read_fd_), wk_fd(wk_fd_)
    {


        //memset(read_buf, 0x00, 256);
        //memset(write_buf, 0x00, 256);
    }

    std::deque<int> * work_list;
    Shannon_C_setting setting;
    pthread_mutex_t * work_lock_ptr;
    int init_total_work;
    int process_i;

    int read_fd;
    int wk_fd;

    char read_buf[256];
    char write_buff[256];
    std::string my_id;

    std::string finish_token;

    static void * run_multi_seq_graph(void * seq_graph_work_ptr)
    {
        Seq_graph_works sgw =
            (*static_cast<Seq_graph_works * >(seq_graph_work_ptr));
        int num;
        int comp_i = DEFAULT_COMP_NUM;
        int64_t peak_mem = DEFAULT_COMP_MEM;

        while(1)
        {
            // take the work
            //pthread_mutex_lock(sgw.work_lock_ptr);
            std::string last_comp_mem = std::to_string(peak_mem);

            // requset
            int offset = 0;
            memcpy(sgw.write_buff, sgw.my_id.c_str(), sgw.my_id.size());
            sgw.write_buff[sgw.my_id.size()] = ',';
            offset += sgw.my_id.size() + 1;
            memcpy(sgw.write_buff + offset, last_comp_mem.c_str(), last_comp_mem.size());
            offset += last_comp_mem.size();
            sgw.write_buff[offset] = ' ';

            sgw.write_buff[offset+1] = '\0';
            if ((num = write(sgw.wk_fd, sgw.write_buff, strlen(sgw.write_buff))) < 0)
				perror("child - write");
			//else
			//	printf("child %d - wrote %d bytes %s \n", sgw.process_i,  num, sgw.write_buff);

            // receive
            std::string receive_token;
			do
			{
                memset(sgw.read_buf, 0x00, 256);
				if ((num = read(sgw.read_fd, sgw.read_buf, sizeof(sgw.read_buf))) < 0)
					perror("child  - read");
				else
				{
                    receive_token = sgw.read_buf;
					//printf("child get: %s\n", sgw.read_buf);
					break;
				}
			} while (num > 0);

            if(receive_token == sgw.finish_token)
            {
                //std::cout <<"releasing a thread " << std::endl;
                //pthread_mutex_unlock(sgw.work_lock_ptr);
                return ((void*) 0);
            }

            comp_i = boost::lexical_cast<int>(receive_token);
            //pthread_mutex_unlock(sgw.work_lock_ptr);

            Sequence_graph_handler * seq_graph_ptr = new Sequence_graph_handler(sgw.setting);

            //if(num_comp_left%SEQ_GRAPH_WORK_STEP_UPDATE==0)
            //{
            //    int percentage = (100 * (seq_graph_work.init_total_work-num_comp_left))/seq_graph_work.init_total_work;
            //    printf("[ %3d%%] MB process %4d, \t num left %d/%d\n",
            //            percentage, seq_graph_work.process_i, num_comp_left,
            //            seq_graph_work.init_total_work);
            //}
            //std::cout << "process  " << seq_graph_work.process_i << " runs comp "
            //          << comp_i << std::endl;
            peak_mem = seq_graph_ptr->run_it(comp_i, false);

            delete seq_graph_ptr;

        }

        return NULL;
    }
};
#endif	/* SEQUENCE_GRAPH_HANDLER_H */
