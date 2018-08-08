/*
 * File:   main.cpp
 * Author: bx
 *
 * Created on August 15, 2017, 2:47 PM
 */
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>

#include "run_tasks.h"

void command_line_for_shannon(int argc, char** argv, Shannon_C_setting & setting);
void command_line_for_partition(int argc, char** argv, Shannon_C_setting & setting);
void command_line_for_seq_graphs(int argc, char** argv, Shannon_C_setting & setting,
                        std::string& kmer_path, std::string& s_read_path,
                        std::string& p1_read_path, std::string& p2_read_path);
void command_line_for_sparse_flow(int argc, char** argv,
            Shannon_C_setting & setting,
            std::string & node_path, std::string & edge_path,
            std::string & path_path, std::string & output_path);
void command_line_for_multi_graph_both(int argc, char** argv, Shannon_C_setting & setting);
void command_line_for_multi_graph_mb(int argc, char** argv, Shannon_C_setting & setting);
void command_line_for_multi_graph_sf(int argc, char** argv, Shannon_C_setting & setting);
void command_line_for_find_rep(int argc, char** argv, Shannon_C_setting & setting,
                        std::string & input_path, std::string & output_path);
void command_line_for_ref_align(int argc, char** argv, Shannon_C_setting & setting,
                        std::string & infile, std::string & ref_file,
                        std::string & outdir);

void add_partition_options(boost::program_options::options_description & desc);
void add_kmer_strand_options(boost::program_options::options_description & desc);
void add_sparse_flow_options(boost::program_options::options_description & desc);
void add_output_options(boost::program_options::options_description & desc);
void add_read_length_options(boost::program_options::options_description & desc);
void add_read_path_options(boost::program_options::options_description & desc);
void add_multi_bridge_options(boost::program_options::options_description & desc);
void add_multi_graph_options(boost::program_options::options_description & desc);

void parse_kmer_strand_options(
    boost::program_options::variables_map & vm, Shannon_C_setting  & setting);
void parse_read_len_and_path(
    boost::program_options::variables_map & vm, Shannon_C_setting  & setting);
void parse_read_len_only(boost::program_options::variables_map & vm, Shannon_C_setting  & setting);
void parse_required_multi_graph(
    boost::program_options::variables_map & vm, Shannon_C_setting  & setting);
void parse_partition_option(
        boost::program_options::variables_map & vm, Shannon_C_setting  & setting);
void parse_seq_graph_optional(
        boost::program_options::variables_map & vm, Shannon_C_setting  & setting);
void parse_sparse_flow_optional(
        boost::program_options::variables_map & vm, Shannon_C_setting  & setting);
void parse_output_options(
        boost::program_options::variables_map & vm, Shannon_C_setting  & setting);

void print_partition_setting(Shannon_C_setting  & setting);

int main(int argc, char** argv) {
    bool show_command_help_message = false;
    if (argc > 1)
    {
        std::string subcommand(argv[1]);
        Shannon_C_setting setting;
        if(subcommand == "custom")
        {
            if(argc != 3)
            {
                std::cout << "\033[1;31m";
                std::cout << "usage: ./Shannon_C_seq custom setting_file_path\n";
                std::cout << "\033[0m";
                exit(1);
            }

            std::string setting_file_name(argv[2]);
            boost::filesystem::path base_path = boost::filesystem::current_path();
            std::string base_path_str(base_path.c_str());
            std::string setting_path = base_path_str + "/" + setting_file_name;

            if(!boost::filesystem::exists(setting_path))
            {
                std::cout << "\033[1;31m";
                std::cout << "usage: ./Shannon_C_seq custom setting_file_name\n";
                std::cout << "but setting_file_path not exist: " << setting_path << std::endl;;
                std::cout << "\033[0m";
                exit(1);
            }
            run_custom(argc, argv, setting_path);
        }
        else if(subcommand == "shannon")
        {
            command_line_for_shannon(argc, argv, setting);
            test_all(&setting);
        }
        else if(subcommand == "partition")
        {
            command_line_for_partition(argc, argv, setting);
            test_Contig_graph(&setting);
        }
        else if(subcommand == "multi-graph")
        {
            bool show_mg_help_message = false;
            if (argc > 2)
            {
                std::string sub_sub_cmd = argv[2];
                if (sub_sub_cmd == "both")
                {
                    command_line_for_multi_graph_both(argc, argv, setting);

                    Multi_graph_handler multi_graph(setting);
                    multi_graph.run_multi_seq_graph();
                    // sparse flow
                    multi_graph.run_multi_sparse_flow(-1);
                }
                else if (sub_sub_cmd == "multi-bridge")
                {
                    command_line_for_multi_graph_mb(argc, argv, setting);
                    test_multi_seq_graph(&setting);
                    //run multi_bridge
                }
                else if (sub_sub_cmd == "sparse-flow")
                {
                    command_line_for_multi_graph_sf(argc, argv, setting);

                    Multi_graph_handler multi_graph(setting);
                    multi_graph.run_multi_sparse_flow(-1);
                }
                else
                {
                    std::cout << "\033[1;31m";
                    std::cout << "unknown mode: " << sub_sub_cmd << std::endl;
                    show_mg_help_message = true;
                    std::cout << "\033[0m";
                }
            }
            else
            {
                show_mg_help_message = true;
            }


            if( show_mg_help_message)
            {
                std::cout << "\033[0;33m";
                std::cout << "after subcommand specify operating mode"  << std::endl;
                std::cout << "  multi-bridge : batch multi-bridge"  << std::endl;
                std::cout << "  sparse-flow  : batch sparse flow"  << std::endl;
                std::cout << "  both         : run multi-bridge first then sparse-flow" << std::endl;
                std::cout << "usage: [subcommand] [mode] [options]" << std::endl;
                std::cout << "usage: [subcommand] [mode] --help to see options for each mode" << std::endl;
                std::cout << "\033[0m";
                exit(0);
            }
        }
        else if(subcommand == "multi-bridge")
        {
            std::string kmer_path;
            std::string s_read_path;
            std::string p1_read_path;
            std::string p2_read_path;
            command_line_for_seq_graphs(argc, argv, setting,
                    kmer_path, s_read_path, p1_read_path, p2_read_path);
            Sequence_graph_handler seq_graph_handler(setting);

            int comp_i = -1;
            seq_graph_handler.setup_input_file_path(kmer_path, s_read_path,
                                                    p1_read_path, p2_read_path);
            seq_graph_handler.run_it(comp_i, true);
        }
        else if(subcommand == "sparse-flow")
        {
            std::string node_path, edge_path, path_path, output_path;
            command_line_for_sparse_flow(argc, argv, setting,
                        node_path, edge_path, path_path, output_path);
            Comp_graph comp_graph(-1,-1);
            pthread_mutex_t write_lock;
            if (pthread_mutex_init(&write_lock, NULL) != 0)
                { shc_log_error("unable to initialize mutex\n"); exit(1);}

            Sparse_flow_handler sparse_flow_handler(setting);
            sparse_flow_handler.run_sparse_flow_handler_helper( comp_graph,
                node_path, edge_path, path_path,
                &write_lock, output_path);
        }
        else if(subcommand == "find-rep")
        {
            std::string infile, outfile;
            command_line_for_find_rep(argc, argv, setting, infile, outfile);
            Multi_graph_handler multi_graph(setting);

            multi_graph.find_representatives(infile, outfile);
        }
        else if(subcommand == "ref-align")
        {
            std::string infile, outdir, ref_file;
            command_line_for_ref_align(argc, argv, setting, infile, ref_file, outdir);
            Local_files & lf = setting.local_files;
            lf.output_path = outdir;
            lf.summary_file_path = outdir + "/summary.log";
            lf.eval_dir_path = outdir;
            lf.eval_path = outdir + "/eval.log";
            lf.reconstructed_seq_path = infile;
            lf.reference_seq_path = ref_file;
            eval_reconstructed_seq(lf);
        }
        else
        {
            std::cout << "\033[1;31m";
            std::cout << "unknown subcommand: " << subcommand << std::endl;
            show_command_help_message = true;
            std::cout << "\033[0m";
        }

        std::cout << "Shannon " << subcommand << " finish" << std::endl << std::endl;
    }
    else
    {
        show_command_help_message = true;
    }

    if(show_command_help_message)
    {
        std::cout << "\033[0;33m";
        std::cout << "\nAvailable subcommands are : " << std::endl;
        std::cout << "custom      : " << " Use config file to provide detailed procedure "<< std::endl;
        std::cout << "shannon     : " << " Run entire Shannon "<< std::endl;
        std::cout << "partition   : " << " Given a large noisy fasta-reads, error correction and partition " << std::endl
                  << "              " << " the reads into components" << std::endl;
        std::cout << "multi-graph : " << " Given a group of fasta-reads, kmers pairs, it allows " << std::endl
                  << "              " << " batch processes on all components. " << std::endl
                  << "              " << " It operates on three modes" << std::endl;
        std::cout << "multi-bridge: " << " Given one fasta-reads, kmer pair, it output multi-bridged " << std::endl
                  << "              " << " de Bruijn graph and singleton nodes" << std::endl;
        std::cout << "sparse-flow : " << " Given a graph stored in the format from multi-bridge, " << std::endl
                  << "              " << " it run sparse flow algorithm to get final output seq" << std::endl;
        std::cout << "find-rep    : " << " Reduce redundant transcripts" << std::endl;
        std::cout << "ref-algin   : " << " Align to reference file, and generate evaulation files" << std::endl;
        std::cout << std::endl;
        std::cout << "usage: [subcommand] [options]" << std::endl;
        std::cout << "usage: [subcommand] --help to see options" << std::endl;
        std::cout << "\033[0m";
        exit(0);
    }

    return 0;
}

void parse_sparse_flow_optional(
        boost::program_options::variables_map & vm, Shannon_C_setting  & setting)
{
    setting.kmer_length = vm["kmer_length"].as<int>();
    setting.num_parallel = vm["num_process"].as<int>();
    setting.sparse_flow_setup.multiple_test = vm["multiple_test"].as<int>();
    setting.output_seq_min_len = vm["output_seq_min_len"].as<int>();
}

void parse_required_multi_graph(
    boost::program_options::variables_map & vm, Shannon_C_setting  & setting)
{
    Local_files & lf = setting.local_files;

    lf.output_path = vm["output_dir"].as<std::string >();
    if(!is_abs_path(lf.output_path))
        convert_relative_path_to_abs(lf.output_path);
    //std::cout << "lf.output_path " << lf.output_path<< std::endl;
    add_directory_if_not_exist(lf.output_path);
    lf.reset_paths();


    if (vm.count("single_read_length"))
    {
        setting.has_single = true;
        setting.single_read_length = vm["single_read_length"].as<int>();
    }

    if (vm.count("pair_read_length"))
    {
        setting.has_pair = true;
        std::vector<int>  pair_lenght =
                vm["pair_read_length"].as<std::vector<int> >();
        if(pair_lenght.size()!=2)
        {
            std::cerr << "pair read length error " << std::endl;
        }
        setting.pair_1_read_length = pair_lenght.at(0);
        setting.pair_2_read_length = pair_lenght.at(1);
    }

    if(!setting.has_pair && !setting.has_single)
    {
        std::cout << "\033[0;33m";
        std::cout << "ERROR: input length needed"  << std::endl;
        std::cout << "\033[0m";
        exit(0);
    }

    lf.output_components_read_dir = vm["read_components_dir"].as<std::string>();
    if(!convert_to_abs_and_check_exist(lf.output_components_read_dir))
        exit(0);

    lf.output_components_kmer_dir = vm["kmer_components_dir"].as<std::string>();
    if(!convert_to_abs_and_check_exist(lf.output_components_kmer_dir))
        exit(0);
}

void command_line_for_multi_graph_mb(int argc, char** argv, Shannon_C_setting & setting)
{
    //since input args are the same
    command_line_for_multi_graph_both(argc, argv, setting);
}

void command_line_for_find_rep(int argc, char** argv, Shannon_C_setting & setting,
                        std::string & input_path, std::string & output_path)
{
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    try {
        desc.add_options()
            ("help", "produce help message")
            ("output_seq_min_len,e", po::value<int>()->default_value(200),
                        "minimal reconstructed output length")
            ("input", po::value<std::string>()->required(),
                        "input path")
            ("output", po::value<std::string>()->required(),
                        "output path")
            ("double_strand,d", po::value<bool>()->default_value(true),
                  "Specify if reads are considered double stranded.")
        ;

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);

        if (vm.count("help")) {
            std::cout << "Usage: options_description [options]\n";
            std::cout << desc;
            exit(0);
        }

        po::notify(vm);
        setting.output_seq_min_len = vm["output_seq_min_len"].as<int>();
        setting.is_double_stranded = vm["double_strand"].as<bool>();
        input_path = vm["input"].as<std::string>();
        if(!convert_to_abs_and_check_exist(input_path))
            exit(0);

        output_path = vm["output"].as<std::string>();
        if(!is_abs_path(output_path))
            convert_relative_path_to_abs(output_path);
        if(exist_path(output_path))
            overwirte_a_file(output_path);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }
}

void command_line_for_multi_graph_sf(int argc, char** argv, Shannon_C_setting & setting)
{
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    try {
        desc.add_options()
            ("help", "produce help message")
            ("graphs_dir,h", po::value<std::string>()->required(),
                        "specify graphs dir generated from multi-bridge")
        ;
        add_output_options(desc);
        add_sparse_flow_options(desc);
        add_multi_graph_options(desc);
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);

        if (vm.count("help")) {
            std::cout << "Usage: options_description [options]\n";
            std::cout << desc;
            exit(0);
        }

        po::notify(vm);

        Local_files & lf = setting.local_files;

        lf.output_path = vm["output_dir"].as<std::string >();
        if(!is_abs_path(lf.output_path))
            convert_relative_path_to_abs(lf.output_path);
        add_directory_if_not_exist(lf.output_path);
        lf.reset_paths();
        setting.output_seq_min_len = vm["output_seq_min_len"].as<int>();

        lf.output_seq_graph_path = vm["graphs_dir"].as<std::string>();
        if(!is_abs_path(lf.output_seq_graph_path))
            convert_relative_path_to_abs(lf.output_seq_graph_path);

        setting.num_parallel = vm["num_process"].as<int>();
        setting.sparse_flow_setup.multiple_test = vm["multiple_test"].as<int>();



    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }

}

void command_line_for_multi_graph_both(int argc, char** argv, Shannon_C_setting & setting)
{
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    try {
        desc.add_options()
            ("help", "produce help message")
            ("kmer_length,k", po::value<int>()->default_value(25),
                  "kmer length. default 25")
            ("double_strand,d", po::value<bool>()->default_value(true),
                  "is reads considered double stranded.")
            ("read_components_dir,r", po::value<std::string>()->required(),
                        "specify reads_components dir produced from partition step")
            ("kmer_components_dir,c", po::value<std::string>()->required(),
                        "specify kmer_components dir produced from partition step")
        ;
        add_read_length_options(desc);
        add_output_options(desc);
        add_sparse_flow_options(desc);
        add_multi_bridge_options(desc);
        add_multi_graph_options(desc);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);

        if (vm.count("help"))
        {
            std::cout << "Usage: options_description [options]\n";
            std::cout << desc;
            exit(0);
        }
        po::notify(vm);

        parse_required_multi_graph(vm, setting);
        setting.kmer_length = vm["kmer_length"].as<int>();
        setting.is_double_stranded = vm["double_strand"].as<bool>();
        setting.num_parallel = vm["num_process"].as<int>();
        setting.seq_graph_setup.max_hop_path = vm["max_hop_for_known_path"].as<int>();
        setting.sparse_flow_setup.multiple_test = vm["multiple_test"].as<int>();
        setting.output_seq_min_len = vm["output_seq_min_len"].as<int>();
        setting.take_single_node_seq = true;
        setting.take_contig_seq = false;

        Local_files & lf = setting.local_files;
        print_and_log_kmer_strand_setting(setting);
        print_and_log_read_length_setting(setting);
        std::cout << "\033[1;33m";
        std::cout << "components_kmer_dir      : " << lf.output_components_kmer_dir << std::endl;
        std::cout << "components_read_dir      : " << lf.output_components_read_dir << std::endl;
        std::cout << "\033[0m";

        print_and_log_output_setting(setting);

        shc_log_info(shc_logname, "output_components_kmer_dir: %s\n", lf.output_components_kmer_dir.c_str());
        shc_log_info(shc_logname, "output_components_read_dir: %s\n", lf.output_components_read_dir.c_str());
        print_and_log_multi_graph_setting(setting);
        print_and_log_mb_setting(setting);
        print_and_log_sf_setting(setting);

    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }
}

void command_line_for_sparse_flow(
        int argc, char** argv, Shannon_C_setting & setting,
        std::string & node_path, std::string & edge_path,
        std::string & path_path, std::string & output_path)
{
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    try {
        desc.add_options()
            ("help", "produce help message")
            ("node_path", po::value<std::string>()->required(),
                                        "input nodes path to sparse flow")
            ("edge_path", po::value<std::string>()->required(),
                                        "input edges path to sparse flow")
            ("path_path", po::value<std::string>()->required(),
                                        "input path path to sparse flow")
            ("output_seq_min_len,e", po::value<int>()->default_value(200),
                      "minimal reconstructed output length")
            ("output_path", po::value<std::string>()->required(),
                                    "output path")
        ;

        add_sparse_flow_options(desc);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);
        if (vm.count("help")) {
            std::cout << "Usage: options_description [options]\n";
            std::cout << desc;
            exit(0);
        }
        po::notify(vm);

        setting.output_seq_min_len = vm["output_seq_min_len"].as<int>();
        setting.sparse_flow_setup.multiple_test = vm["multiple_test"].as<int>();
        output_path = vm["output_path"].as<std::string>();
        if(!is_abs_path(output_path))
            convert_relative_path_to_abs(output_path);
        if(exist_path(output_path))
            overwirte_a_file(output_path);

        node_path = vm["node_path"].as<std::string>();
        edge_path = vm["edge_path"].as<std::string>();
        path_path = vm["path_path"].as<std::string>();
        if (!convert_to_abs_and_check_exist(node_path) ||
            !convert_to_abs_and_check_exist(edge_path) ||
            !convert_to_abs_and_check_exist(path_path) )
        {
            exit(0);
        }
        std::cout << "\033[1;33m";
        std::cout << "results saves to          : " << output_path << std::endl;
        std::cout << "output_seq_min_len        : " << setting.output_seq_min_len << std::endl;
        std::cout << "node_path                 : " << node_path << std::endl;
        std::cout << "edge_path                 : " << edge_path << std::endl;
        std::cout << "path_path                 : " << path_path << std::endl;
        std::cout << "\033[0m";
        shc_log_info(shc_logname, "results saves to          : %s\n", output_path.c_str());
        shc_log_info(shc_logname, "output_seq_min_len        : %d\n", setting.output_seq_min_len);
        shc_log_info(shc_logname, "node_path                 : %s\n", node_path.c_str());
        shc_log_info(shc_logname, "edge_path                 : %s\n", edge_path.c_str());
        shc_log_info(shc_logname, "path_path                 : %s\n", path_path.c_str());
        print_and_log_sf_setting(setting);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }
}

void command_line_for_seq_graphs(int argc, char** argv, Shannon_C_setting & setting,
        std::string& kmer_path, std::string& s_read_path,
        std::string& p1_read_path, std::string& p2_read_path)
{
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    try {
        desc.add_options()
            ("help", "produce help message")
            ("kmer_path,f", po::value<std::string>()->required(),
                      "kmer path used as trusted debruijn graph edge")
        ;

        add_read_length_options(desc);
        add_read_path_options(desc);
        add_multi_bridge_options(desc);
        add_output_options(desc);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);
        if (vm.count("help")) {
            std::cout << "Usage: options_description [options]\n";
            std::cout << desc;
            exit(0);
        }

        po::notify(vm);

        setting.seq_graph_setup.max_hop_path =
                                        vm["max_hop_for_known_path"].as<int>();

        parse_output_options(vm, setting);
        kmer_path = vm["kmer_path"].as<std::string>();
        if(!convert_to_abs_and_check_exist(kmer_path))
            exit(0);
        setting.kmer_length = get_kmer_length_from_kmer_file(kmer_path);

        parse_read_len_and_path(vm, setting);
        s_read_path = setting.local_files.input_read_path;
        p1_read_path = setting.local_files.input_read_path_1;
        p2_read_path = setting.local_files.input_read_path_2;


        Local_files & lf = setting.local_files;
        print_and_log_kmer_strand_setting(setting);
        print_and_log_read_length_setting(setting);
        std::cout << "\033[1;33m";
        std::cout << "kmer_path:                " << kmer_path << std::endl;
        std::cout << "\033[0m";
        print_and_log_read_length_setting(setting);
        print_and_log_input_path_setting(setting);
        print_and_log_output_setting(setting);
        print_and_log_mb_setting(setting);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }

}


/*
Always parse output first, the output file tree is then setup by using reset
command, then input path should be set to override the auto-generated path */
void parse_output_options(
        boost::program_options::variables_map & vm, Shannon_C_setting  & setting)
{
    setting.output_seq_min_len = vm["output_seq_min_len"].as<int>();
    setting.local_files.output_path = vm["output_dir"].as<std::string>();

    if(!is_abs_path(setting.local_files.output_path))
    {
        convert_relative_path_to_abs(setting.local_files.output_path);
    }

    add_directory_if_not_exist(setting.local_files.output_path);
    setting.local_files.reset_paths();
}

void parse_seq_graph_optional(
        boost::program_options::variables_map & vm, Shannon_C_setting  & setting)
{
    setting.kmer_length = vm["kmer_length"].as<int>();
    setting.num_parallel = vm["num_process"].as<int>();
    setting.seq_graph_setup.max_hop_path = vm["max_hop_for_known_path"].as<int>();
}

void add_sparse_flow_options(boost::program_options::options_description & desc)
{
    namespace po = boost::program_options;
    desc.add_options()
        ("multiple_test", po::value<int>()->default_value(8),
                "number of times that a linear program sparse flow is solved")
    ;
}

void add_output_options(boost::program_options::options_description & desc)
{
    namespace po = boost::program_options;
    desc.add_options()
        ("output_seq_min_len,e", po::value<int>()->default_value(200),
                  "minimal reconstructed output length")
        ("output_dir,o", po::value<std::string>()->required(),
                              "output dir path")
    ;
}

void add_read_length_options(boost::program_options::options_description & desc)
{
    namespace po = boost::program_options;
    desc.add_options()
        ("single_read_length,l", po::value<int>(),
                                          "single read length")
        ("pair_read_length,i", po::value<std::vector<int> >()->multitoken(),
                                          "pair read length")
    ;
}

void add_read_path_options(boost::program_options::options_description & desc)
{
    namespace po = boost::program_options;
    desc.add_options()
        ("SE_read_path,s", po::value<std::string>(),
                    "single ended read path ")
        ("PE_read_path,p", po::value<std::vector<std::string> >()->multitoken(),
                                          "pair read path")
    ;
}

void add_multi_bridge_options(boost::program_options::options_description & desc)
{
    namespace po = boost::program_options;
    desc.add_options()
        ("max_hop_for_known_path", po::value<int>()->default_value(30),
                "Max number of hops allowable for a read to cover")
    ;
}

void add_multi_graph_options(boost::program_options::options_description & desc)
{
    namespace po = boost::program_options;
    desc.add_options()
        ("num_process,t",po::value<int>()->default_value(1), "Specify number processs for multi graph procedures"
                     "where each process produces two threads by default")
    ;
}

void parse_kmer_strand_options(
    boost::program_options::variables_map & vm, Shannon_C_setting  & setting)
{
    parse_read_len_and_path(vm, setting);
    setting.kmer_length = vm["kmer_length"].as<int>();
    setting.is_double_stranded = vm["double_strand"].as<bool>();
}

void add_kmer_strand_options(boost::program_options::options_description & desc)
{
    namespace po = boost::program_options;
    desc.add_options()
        ("kmer_length,k", po::value<int>()->default_value(25),
              "kmer length. default 25")
        ("double_strand,d", po::value<bool>()->default_value(true),
              "is reads considered double stranded.")
    ;
}

void add_partition_options(boost::program_options::options_description & desc)
{
    namespace po = boost::program_options;
    desc.add_options()
        ("jf_path", po::value<std::string>(),
                  "jf path from jellyfish")
        ("kmer_path", po::value<std::string>(),
                    "kmer path from jellyfish")
        ("load_factor", po::value<double>()->default_value(0.8),
                                "kmer dictionary load factor, control memory usage")
        ("num_sort_thread,u",po::value<int>()->default_value(1),
                        "argument passing to linux sort function")
        ("rmer_length", po::value<int>()->default_value(15),
                        "rmer length for error correction")
        ("threshold", po::value<double>()->default_value(0.5),
                        "ratio of common rmer between two contigs for determing "
                        "if two contifs are repetative")
        ("min_count", po::value<int>()->default_value(3),
                        "minimal kmer counts requirement for serving as seeds "
                        "in the greedy contig formation step, also used in "
                        "further selection formula discussed in the paper")
        ("min_length", po::value<int>()->default_value(75),
                        "minimal contigs length, also also used in "
                        "further selection formula discussed in the paper")
        ("is_use_set", po::value<bool>()->default_value(false),
                        "use accummulating set for error correction, "
                        "it would reduce memory usage, but give less performance")
        ("read_num_test", po::value<int>()->default_value(3),
                        "number test for a read to decide which components "
                        "the read should go to")
        ("is_assign_best",po::value<bool>()->default_value(true),
                        "only assign a read to its best matching component "
                        "or all components which work")
        ("read_sampler_k,g", po::value<int>()->default_value(0),
                        "given that many reads are redundant, setting k to "
                        "subsample matching reads to speedup later works, "
                        "where k can be interpreted as coverage depth. "
                        "k=0: no subsample; k>0: probabilistic sample; "
                        "k<0: sample every (-k) read")
        ("metis_use_multiple_partition",po::value<bool>()->default_value(true),
                        "after one metis partition, reweight graph so that "
                        "previously broken edge stays, it makes sure no "
                        "information loss due to partition, but generate more output")
        ("metis_partition_size",po::value<int>()->default_value(500),
                        "metis parameter to control number of contig within "
                        "partitioned components")
        ("metis_non_partition_size",po::value<int>()->default_value(500),
                        "for those graphs small enough not to be partitioned "
                        "Shannon groups little pieces together, and the parameter "
                        "control the collector component size")
        ("penalty", po::value<int>()->default_value(5),
                        "metis parameter, for controlling partition methods")
        ("overload", po::value<double>()->default_value(2.0),
                        "metis parameter, for controlling partition methods")
    ;
}

void command_line_for_shannon(int argc, char** argv, Shannon_C_setting  & setting)
{
    namespace po = boost::program_options;
    Local_files & lf = setting.local_files;
    try {
        std::string config_path;

        po::options_description desc("Allowed options");
        desc.add_options() ("help", "produce help message");
        add_kmer_strand_options(desc);
        add_read_length_options(desc);
        add_read_path_options(desc);
        add_output_options(desc);
        add_partition_options(desc);

        desc.add_options()
            ("reference,r", po::value<std::string>(),
                                     "Reference path for the output to compare with")
        ;
        add_multi_graph_options(desc);
        add_multi_bridge_options(desc);
        add_sparse_flow_options(desc);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);


        if (vm.count("help")) {
            std::cout << "Usage: options_description [options]\n";
            std::cout << desc;
            exit(0);
        }
        po::notify(vm);

        parse_kmer_strand_options(vm, setting);
        parse_output_options(vm, setting);
        parse_read_len_and_path(vm, setting);
        parse_partition_option(vm, setting);

        setting.num_parallel = vm["num_process"].as<int>();
        setting.seq_graph_setup.max_hop_path = vm["max_hop_for_known_path"].as<int>();
        setting.sparse_flow_setup.multiple_test = vm["multiple_test"].as<int>();
        if(vm.count("reference"))
            lf.reference_seq_path =  vm["reference"].as<std::string>();

        print_and_log_kmer_strand_setting(setting);
        print_and_log_read_length_setting(setting);
        print_and_log_input_path_setting(setting);
        print_and_log_output_setting(setting);
        print_and_log_partition_setting(setting);
        print_and_log_multi_graph_setting(setting);
        print_and_log_mb_setting(setting);
        print_and_log_sf_setting(setting);

        return;
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }
}

void command_line_for_partition(int argc, char** argv, Shannon_C_setting & setting)
{
    namespace po = boost::program_options;
    Local_files & lf = setting.local_files;
    try {
        po::options_description desc("Allowed options");
        desc.add_options() ("help", "produce help message");
        add_kmer_strand_options(desc);
        add_read_length_options(desc);
        add_read_path_options(desc);
        add_output_options(desc);
        add_partition_options(desc);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);

        if (vm.count("help")) {
            std::cout << "Usage: options_description [options]\n";
            std::cout << desc;
            exit(0);
        }

        po::notify(vm);

        parse_read_len_and_path(vm, setting);
        parse_kmer_strand_options(vm, setting);
        parse_output_options(vm, setting);
        parse_partition_option(vm, setting);

        print_and_log_kmer_strand_setting(setting);
        print_and_log_read_length_setting(setting);
        print_and_log_input_path_setting(setting);
        print_and_log_output_setting(setting);
        print_and_log_partition_setting(setting);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }
}

void parse_partition_option(
        boost::program_options::variables_map & vm, Shannon_C_setting  & setting)
{
    Local_files & lf = setting.local_files;
    bool is_provided_jf = false;
    if (vm.count("jf_path"))
    {
        is_provided_jf = true;
        lf.input_jf_path = vm["jf_path"].as<std::string>();
        if(!is_abs_path(lf.input_jf_path))
            convert_relative_path_to_abs(lf.input_jf_path);
    }
    if (vm.count("kmer_path"))
    {
        if(!is_provided_jf)
        {
            std::cout << "\033[0;33m";
            std::cout << "ERROR: jf path must be provided with kmer path" << std::endl;
            std::cout << "\033[0m";
            exit(0);
        }
        lf.input_kmer_path = vm["kmer_path"].as<std::string>();
        if(!convert_to_abs_and_check_exist(lf.input_kmer_path))
            exit(0);
    }

    setting.kmer_length = vm["kmer_length"].as<int>();
    setting.is_double_stranded = vm["double_strand"].as<bool>();
    setting.output_seq_min_len = vm["output_seq_min_len"].as<int>();
    // dup setting
    setting.dup_setting.num_sort_thread = vm["num_sort_thread"].as<int>();

    setting.dup_setting.rmer_length = vm["rmer_length"].as<int>();

    setting.dup_setting.threshold = vm["threshold"].as<double>();

    setting.dup_setting.min_count = vm["min_count"].as<int>();

    setting.dup_setting.min_len = vm["min_length"].as<int>();

    setting.dup_setting.is_use_set = vm["is_use_set"].as<bool>();

    setting.dup_setting.load_factor = vm["load_factor"].as<double>();
// contig graph dup_setting

    setting.contig_graph_setup.num_test = vm["read_num_test"].as<int>();

    setting.contig_graph_setup.is_assign_best = vm["is_assign_best"].as<bool>();

    setting.contig_graph_setup.read_sampler_k = vm["read_sampler_k"].as<int>();
// metis metis_setup

    setting.metis_setup.is_multiple_partition = vm["metis_use_multiple_partition"].as<bool>();

    setting.metis_setup.partition_size = vm["metis_partition_size"].as<int>();

    setting.metis_setup.non_partition_size = vm["metis_non_partition_size"].as<int>();

    setting.metis_setup.penalty = vm["penalty"].as<int>();

    setting.metis_setup.overload = vm["overload"].as<double>();
}

void parse_read_len_only(boost::program_options::variables_map & vm, Shannon_C_setting  & setting)
{
    bool has_se_len = false;
    if (vm.count("single_read_length"))
    {
        setting.has_single = true;
        setting.single_read_length = vm["single_read_length"].as<int>();
        setting.local_files.has_single = true;
        has_se_len = true;
    }

    bool has_pe_len = false;
    if (vm.count("pair_read_length"))
    {
        setting.has_pair = true;
        setting.local_files.has_pair = true;
        std::vector<int>  pair_lenght =
                vm["pair_read_length"].as<std::vector<int> >();
        if(pair_lenght.size()!=2)
        {
            std::cerr << "pair read length error " << std::endl;
        }
        setting.pair_1_read_length = pair_lenght.at(0);
        setting.pair_2_read_length = pair_lenght.at(1);
        has_pe_len = true;
    }

    if(!has_se_len && !has_pe_len)
    {
        std::cout << "\033[0;33m";
        std::cout << "ERROR: input length must be provided"  << std::endl;
        std::cout << "\033[0m";
        exit(0);
    }
}

void parse_read_len_and_path(
    boost::program_options::variables_map & vm, Shannon_C_setting  & setting)
{
    Local_files & lf = setting.local_files;
    bool se_len = false;
    bool se_path = false;
    bool pe_len = false;
    bool pe_path = false;

    if (vm.count("single_read_length"))
    {
        setting.has_single = true;
        setting.single_read_length = vm["single_read_length"].as<int>();
        setting.local_files.has_single = true;
        se_len = true;
    }

    if (vm.count("pair_read_length"))
    {
        setting.has_pair = true;
        setting.local_files.has_pair = true;
        std::vector<int>  pair_lenght =
                vm["pair_read_length"].as<std::vector<int> >();
        if(pair_lenght.size()!=2)
        {
            std::cerr << "pair read length error " << std::endl;
        }
        setting.pair_1_read_length = pair_lenght.at(0);
        setting.pair_2_read_length = pair_lenght.at(1);
        pe_len = true;
    }

    if(!se_len && !pe_len)
    {
        std::cout << "\033[0;33m";
        std::cout << "ERROR: input length must be provided"  << std::endl;
        std::cout << "\033[0m";
        exit(0);
    }

    if (vm.count("SE_read_path"))
    {
        lf.input_read_path = vm["SE_read_path"].as< std::string >();
        if(!convert_to_abs_and_check_exist(lf.input_read_path))
            exit(0);

        lf.has_single = true;
        se_path = true;
    }

    if (vm.count("PE_read_path")) {
        lf.has_pair = true;
        lf.input_read_path_1 = vm["PE_read_path"].as< std::vector<std::string> >().at(0);
        lf.input_read_path_2 = vm["PE_read_path"].as< std::vector<std::string> >().at(1);

        if(!convert_to_abs_and_check_exist(lf.input_read_path_1) ||
           !convert_to_abs_and_check_exist(lf.input_read_path_2) )
        {
            exit(0);
        }

        pe_path = true;
    }

    if(!pe_path && !se_path)
    {
        std::cout << "\033[0;33m";
        std::cout << "ERROR: input path must be provided"  << std::endl;
        std::cout << "\033[0m";
        exit(0);
    }

    if(pe_path && pe_len || !pe_path && !pe_len)
        ;
    else
    {
        std::cout << "\033[0;33m";
        std::cout << "ERROR: pair length and pair path must be provided at the same time" << std::endl;
        std::cout << "\033[0m";
        exit(0);
    }

    if(se_path && se_len || !se_path && !se_len)
        ;
    else
    {
        std::cout << "\033[0;33m";
        std::cout << "ERROR: single length and single path must be provided at the same time" << std::endl;
        std::cout << "\033[0m";
        exit(0);
    }
}

void command_line_for_ref_align(int argc, char** argv, Shannon_C_setting & setting,
                        std::string & input_path, std::string & ref_file,
                        std::string & output_dir)
{
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    try {
        desc.add_options()
            ("help", "produce help message")
            ("input", po::value<std::string>()->required(),
                        "reconstucted fasta path")
            ("ref", po::value<std::string>()->required(),
                        "reference path")
            ("output_dir", po::value<std::string>()->required(),
                        "output dir")
        ;

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);

        if (vm.count("help")) {
            std::cout << "Usage: options_description [options]\n";
            std::cout << desc;
            exit(0);
        }

        po::notify(vm);
        input_path = vm["input"].as<std::string>();
        if(!convert_to_abs_and_check_exist(input_path))
            exit(0);


        output_dir = vm["output_dir"].as<std::string>();
        if(!convert_to_abs_and_check_exist(output_dir))
            add_directory_if_not_exist(output_dir);

        ref_file = vm["ref"].as<std::string>();
        if(!convert_to_abs_and_check_exist(ref_file))
            exit(0);
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }
}
