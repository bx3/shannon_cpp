#ifndef MULTI_GRAPH_HANDLER_H
#define MULTI_GRAPH_HANDLER_H

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <vector>
#include <queue>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <random>

#include "Sequence_graph_handler.h"
#include "Sparse_flow_handler.h"
#include "shc_type.h"
#include "json_parser.h"
#include "log.h"
#include "local_file_structure.h"

#define SHOW_STEP 10000
#define VEC_INIT_SIZE 10

#define ROUGH_EMPIRICAL_KMER_NUM_MEM_RATIO 600


class Multi_graph_handler {
    friend class Contig_handler;
public:
    typedef std::pair<int, int> ID_size_pair;
    struct Comp_size_info {
        int comp_id;
        int64_t kmer_num;
        int64_t read_size;
        int64_t mem_required;

        Comp_size_info() {}
        Comp_size_info(int comp_id_, int64_t kmer_num_, int64_t read_size_) :
            comp_id(comp_id_), kmer_num(kmer_num_), read_size(read_size_) {}


        bool operator<(const Comp_size_info& p2) const {
    	       return( kmer_num*ROUGH_EMPIRICAL_KMER_NUM_MEM_RATIO + read_size
                        > p2.kmer_num*ROUGH_EMPIRICAL_KMER_NUM_MEM_RATIO + p2.read_size);
        }
    };

    struct Process_comp_man {
        Process_comp_man () {}
        Process_comp_man (int num_parallel_) : num_parallel(num_parallel_),
                        process_comp(num_parallel_)
        {
            kmer_num_mem_ratio_list.push_back(ROUGH_EMPIRICAL_KMER_NUM_MEM_RATIO);
        }

        double get_mean_ratio()
        {
            double sum_ratio = 0.0;
            for(int i=0; i<kmer_num_mem_ratio_list.size(); i++)
            {
                sum_ratio += kmer_num_mem_ratio_list[i];
            }
            return sum_ratio/static_cast<double>(kmer_num_mem_ratio_list.size());
        }

        int num_parallel;
        std::vector<Comp_size_info > process_comp;
        std::vector<double> kmer_num_mem_ratio_list;
    };

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

    //typedef tsl::hopscotch_map<uint64_t, std::list<Seq_info>,
    //               hash_u64, equ64> Rmer_contig_map;
    //typedef tsl::hopscotch_map<uint64_t, std::list<Seq_info>,
    //                hash_u64, equ64>::iterator Rmer_contig_map_iterator;

    typedef spp::sparse_hash_map<uint64_t, std::vector<Seq_info>,
                   hash_u64, equ64> Rmer_contig_map;
    typedef spp::sparse_hash_map<uint64_t, std::vector<Seq_info>,
                   hash_u64, equ64>::iterator Rmer_contig_map_iterator;
    typedef spp::sparse_hash_map<std::string, std::string> Seqs_header_map;
    typedef spp::sparse_hash_map<std::string, std::string>::iterator Seqs_header_map_iterator;


    typedef tsl::hopscotch_map<uint64_t, std::string,
                   hash_u64, equ64> Header_seq_map;
    typedef tsl::hopscotch_map<uint64_t, std::string,
                hash_u64, equ64>::iterator Header_seq_map_iterator;

    Multi_graph_handler(Shannon_C_setting & setting_): setting(setting_)
    {
        mem_safe_ratio = 1.5;
        avail_mem = setting.avail_mem;
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
    void sort_work_list_by_size(std::list<Comp_size_info> & work_list);
    bool get_read_file_size(int i, int64_t & size);
    bool get_comp_num_kmer(int i, int64_t & size);
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
    uint64_t get_parallel_total_mem(std::vector<pid_t> & child_pids, int num_parallel);
    Comp_size_info get_next_comp(std::vector<pid_t> & child_pids,
                        std::list<Comp_size_info> & work_list,
                        int num_parallel, double mean_ratio);
    uint64_t approx_mem_from_kmer(uint64_t kmer_kb);
    bool is_mem_enough_for_comp(Comp_size_info & comp_info, double mean_ratio,
                                int64_t avail_memory, double memory_safe_ratio);

    pthread_mutex_t work_lock;
    pthread_mutex_t write_lock;

    Rmer_contig_map rmer_contig_map;

    int num_comp;
    Shannon_C_setting & setting;
    Block_timer timer;

    int rmer_length;

    double pair_fp_thresh;

    double mem_safe_ratio;

    int64_t total_allowed_mem;
    int64_t avail_mem;
};

#endif
