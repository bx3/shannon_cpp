#ifndef MULTI_GRAPH_HANDLER_H
#define MULTI_GRAPH_HANDLER_H

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <pthread.h>
#include "Sequence_graph_handler.h"
#include "Sparse_flow_handler.h"
#include "shc_type.h"
#include "json_parser.h"
#include "log.h"
#include "local_file_structure.h"
#include <vector>
#include <list>
#include <deque>
#include <algorithm>
#include <random>

#define SHOW_STEP 10000
#define VEC_INIT_SIZE 10


class Multi_graph_handler {
    friend class Contig_handler;
public:
    typedef std::pair<int, size_t> ID_size_pair;
    struct File_sorter {
        bool operator()(const ID_size_pair& p1, const ID_size_pair& p2) const {
    	return( p1.second > p2.second);
        }
    };

    struct Seq_info {
        Seq_info(){}
        Seq_info(const uint64_t & header_, const int64_t & index_) :
                    header(header_), index(index_) {}
        uint64_t header;
        int64_t index;
    };

    typedef tsl::hopscotch_map<uint64_t, std::vector<Seq_info>,
                   hash_u64, equ64> Rmer_contig_map;
    typedef tsl::hopscotch_map<uint64_t, std::vector<Seq_info>,
                    hash_u64, equ64>::iterator Rmer_contig_map_iterator;

    typedef tsl::hopscotch_map<uint64_t, std::string,
                   hash_u64, equ64> Header_seq_map;
    typedef tsl::hopscotch_map<uint64_t, std::string,
                hash_u64, equ64>::iterator Header_seq_map_iterator;

    Multi_graph_handler(Shannon_C_setting & setting_): setting(setting_)
    {
        pair_fp_thresh = 0.9;
        rmer_length = setting_.rmer_length;
        if (pthread_mutex_init(&work_lock, NULL) != 0)
        {
            shc_log_error("unable to initialize mutex\n");
            exit(1);
        }
        if (pthread_mutex_init(&write_lock, NULL) != 0)
        {
            shc_log_error("unable to initialize mutex\n");
            exit(1);
        }
        if(!setting.local_files.single_node_dir.empty())
        {
            //std::cout <<"single_node_dir " << setting.local_files.single_node_dir << std::endl;
            //boost::filesystem::path dir_boost_path(setting.local_files.single_node_dir);
            add_or_overwrite_directory(setting.local_files.single_node_dir,
                                            setting.local_files.output_path);

        }
    }


    int count_num_component_seq_graph();
    int count_num_component_sparse_flow();

    int get_num_components_seq_graph() {return count_num_component_seq_graph();}
    int get_num_components_sparse_flow() {return count_num_component_sparse_flow();}

    void run_multi_seq_graph();
    void run_multi_sparse_flow(int specific_comp);

    void process_multi_seq_graph(int num_process, int num_components);
    void partition_work_to_process(int num_process, std::deque<int> & work_list,
                                std::vector<std::deque<int>> & process_queue);
    void partition_work_to_process_randomize(int num_process, std::deque<Comp_graph> & work_list,
                                std::vector<std::deque<Comp_graph>> & process_queue);

    void process_sparse_flow_graph(int & num_parallel, int num_components, int specific_comp);
    void collect_process_file(int num_parallel);
    void combine_sf_single_seq_output();
    void combine_contigs_seq_output();

    void filter_output_seq_by_length();

    void combine_all_reconst_seq_to_output();

    bool output_fasta_file_validator(std::string & output_path);

    void dump_all_single_nodes_to_reconstruct_output();

    void sort_work_list_by_size(std::deque<int> & work_list);
    bool get_read_file_size(int i, size_t & size);
    // post processing
    void find_representatives(std::string in_file, std::string output_file);
    bool duplicate_check_ends(std::vector<std::string> & header_seq_map, uint64_t header, bool rc);

    // use pair end read to further filtering
    void filter_FP(std::string reconstructed_seq_path);
    void write_filtered_tr(std::string depth_file, std::string in_tr_file,
                           std::string out_tr_file, std::string log_file);



    void filter_reconstructed();
    void deallocate_mem();

private:
    void reverse_complement(std::string & seq);
    void check_if_required_files_exist();

    pthread_mutex_t work_lock;
    pthread_mutex_t write_lock;

    Rmer_contig_map rmer_contig_map;

    int num_comp;
    Shannon_C_setting & setting;
    Block_timer timer;

    int rmer_length;

    double pair_fp_thresh;


};

#endif
