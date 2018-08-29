/*
 * File:   shannon_C_seq_helper.h
 * Author: bx
 *
 * Created on December 1, 2017, 2:45 PM
 */

#ifndef SHANNON_C_SEQ_HELPER_H
#define	SHANNON_C_SEQ_HELPER_H

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <algorithm>
#include <iterator>

#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sys/wait.h>
#include <sys/types.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include "log.h"
#include "json_parser.h"
#include "local_file_structure.h"
#include "shc_type.h"

uint64_t get_mem(pid_t id);
int64_t get_machine_physical_limit_mem();

struct JF_stats {
    JF_stats() {}
    JF_stats(uint64_t uniq_num_, uint64_t dist_num_, uint64_t total_num_,
            uint64_t max_count_) : uniq_num(uniq_num_), dist_num(dist_num_),
            total_num(total_num_), max_count(max_count_) {}
    uint64_t uniq_num;
    uint64_t dist_num;
    uint64_t total_num;
    uint64_t max_count;
};

void fasta_file_validator(std::string path);
void run_jellyfish(Shannon_C_setting & setting);
void run_jellyfish_for_file(Shannon_C_setting & setting, std::string & kmer_path);
void run_command(std::string cmd, bool print_cmd);
void print_yellow_cmd(std::string cmd);

void parse_jf_info(Local_files & lf, JF_stats & jf_info);

void check_read_type_consistency(Shannon_C_setting & setting);
void modify_setting_with_command(Local_files & lf, boost::program_options::variables_map & vm);
int parse_main_command_line(int ac, char** av, Shannon_C_setting & setting);
uint64_t get_num_read(std::string file);

void eval_reconstructed_seq(Local_files & lf);
void eval_reconstructed_seq(std::string reference_seq_path,
                            std::string reconstructed_seq_path,
                            std::string eval_dir_path,
                            std::string eval_path);

std::string get_py_eval_summary(Local_files & lf, std::vector<std::string> & summaries);
void count_seq_types_num(std::string filename, int & num_single_seq,
                        int & num_contig_seq, int & num_sf_seq);
int eval_align_counts(std::string filename, int & num_single_seq,
                        int & num_contig_seq, int & num_sf_seq,
                        int & num_single_contribute, int & num_contig_contribute,
                        int & num_sf_contribute, struct Local_files & lf);

void print_and_log_all_setting(Shannon_C_setting & setting);
void print_and_log_general_setting(Shannon_C_setting & setting);
void print_and_log_partition_setting(Shannon_C_setting & setting);
void print_and_log_multi_graph_setting(Shannon_C_setting & setting);
void print_and_log_mb_setting(Shannon_C_setting & setting);
void print_and_log_sf_setting(Shannon_C_setting & setting);
void print_and_log_general_setting(Shannon_C_setting & setting);
void print_and_log_partition_setting(Shannon_C_setting & setting);
void print_and_log_general_setting(Shannon_C_setting & setting);
void print_and_log_multi_graph_setting(Shannon_C_setting & setting);
void print_and_log_mb_setting(Shannon_C_setting & setting);
void print_and_log_sf_setting(Shannon_C_setting & setting);
void print_and_log_kmer_strand_setting(Shannon_C_setting & setting);
void print_and_log_input_path_setting(Shannon_C_setting & setting);
void print_and_log_read_length_setting(Shannon_C_setting & setting);
void print_and_log_kmer_strand_setting(Shannon_C_setting & setting);
void print_and_log_output_setting(Shannon_C_setting & setting);
void print_and_log_find_rep_setting(Shannon_C_setting & setting);

pid_t fork_mem_profiler(pid_t running_pid, std::string syrupy_path);
void join_mem_profiler(pid_t profiler_pid);


void get_mem_statistics(int num_parallel, Mem_profiler & mp);
uint64_t get_max_RSS(std::string filename);


int get_kmer_length_from_kmer_file(std::string kmer_path);
void produce_summary_file(struct Local_files & lf);

#ifdef USE_APPROX_COUNT
inline kmer_count_t encode_count(uint64_t count, double cr)
{
    if(count <= EXACT_THRESH || cr==1.0)
        return count;
    else
    {
        uint64_t out = EXACT_THRESH + (count-EXACT_THRESH)/cr;
        //std::cout << "out " << out << std::endl;
        assert(out<=KMER_MAX_COUNT);
        return out;
    }
}

inline uint64_t get_count(uint64_t encoded_count, double cr)
{
    if(encoded_count <= EXACT_THRESH || cr==1.0)
        return encoded_count;
    else
    {
        return EXACT_THRESH + (encoded_count-EXACT_THRESH)*cr;
    }
}
#endif

#ifdef USE_EXACT_COUNT
inline kmer_count_t encode_count(uint64_t count, double cr)
{
    return count;
}

inline uint64_t get_count(uint64_t encoded_count, double cr)
{
    return encoded_count;
}
#endif


// code found at
// https://codereview.stackexchange.com/questions/186535/progress-bar-in-c
// with very little change
class Progress_bar
{
    static const auto overhead = sizeof " [100%]";

    std::ostream& os;
    const std::size_t bar_width;
    std::string message;
    const std::string full_bar;
    struct Block_timer bt;

 public:
    Progress_bar(std::ostream& os, std::size_t line_width,
                 std::string message_, const char symbol = '.')
        : os{os},
          bar_width{line_width - overhead},
          message{std::move(message_)},
          full_bar{std::string(bar_width, symbol) + std::string(bar_width, ' ')}
    {
        if (message.size()+1 >= bar_width || message.find('\n') != message.npos) {
            os << message << '\n';
            message.clear();
        } else {
            message += ' ';
        }
        clock_gettime(CLOCK_MONOTONIC, &bt.nano_start);
        write(0.0);
    }

    // not copyable
    Progress_bar(const Progress_bar&) = delete;
    Progress_bar& operator=(const Progress_bar&) = delete;

    ~Progress_bar()
    {
        write(1.0);
        os << '\n';
    }

    void write(double fraction)
    {
        // clamp fraction to valid range [0,1]
        if (fraction < 0)
        fraction = 0;
        else if (fraction > 1)
        fraction = 1;

        auto width = bar_width - message.size();
        auto offset = bar_width - static_cast<unsigned>(width * fraction);

        clock_gettime(CLOCK_MONOTONIC,&bt.nano_stamp);
        bt.nTime = (bt.nano_stamp.tv_sec - bt.nano_start.tv_sec);

        os << '\r' << message;
        os.write(full_bar.data() + offset, width);
        os << " [" << std::setw(3) << static_cast<int>(100*fraction) << "%] "
            <<bt.nTime << "sec"<< std::flush;
        std::cout.flush();
    }

};


#endif	/* SHANNON_C_SEQ_HELPER_H */
