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

void parse_arg_for_shannon(int argc, char** argv, Shannon_C_setting & setting);

int main(int argc, char** argv) {
    std::string subcommand(argv[1]);

    if(subcommand == "custom")
    {
        run_custom(argc, argv);
    }
    else if(subcommand == "shannon")
    {
        Shannon_C_setting setting;
        parse_arg_for_shannon(argc, argv, setting);
        print_setting(setting);
        log_setting(setting);
        test_all(&setting);
    }
    //else if(subcommand == "partition")
    //{
    //    ;
    //}
    //else if(subcommand == "multi_graph")
    //{
    //    ;
    //}
    //else if(subcommand == "seq_graph")
    //{
    //    ;
    //}
    //else if(subcommand == "sparse_flow")
    //{
    //    ;
    //}
    else
    {
        std::cout << "unknown subcommand: " << subcommand << std::endl;
        std::cout << "Available subcommands are : " << std::endl;
        std::cout << "custom      : " << " provide detailed procedure "<< std::endl;
        std::cout << "shannon     : " << " run entire process "<< std::endl;
        //std::cout << "multi-graph : " << " given many fasta-reads(pair), kmers, output reconstructed seq" << std::endl;
        //std::cout << "seq-graph   : " << " given one fasta-reads(pair), kmer, output reconstructed seq" << std::endl;
        //std::cout << "sparse-flow : " << " given a graph stored in files, output reconstructed seq" << std::endl;
        //std::cout << "partition   : " << " given a large noisy fasta-reads(pair), partition and filter the reads " << std::endl;
        std::cout << std::endl;
        std::cout << "usage: [subcommand] [options]" << std::endl;
        exit(0);
    }
    return 0;
}

void parse_arg_for_shannon(int argc, char** argv, Shannon_C_setting  & setting)
{
    namespace po = boost::program_options;
    Local_files & lf = setting.local_files;
    try {
        std::string config_path;

        bool se_len = false;
        bool se_path = false;
        bool pe_len = false;
        bool pe_path = false;

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("shannon", po::value<std::string>(),
                  "")
            ("config,j", po::value<std::string>(),
                  "json config file as default, (required)")
            ("kmer_length,k", po::value<int>()->default_value(25),
                  "specify kmer length.")
            ("single_strand,d", po::value<bool>(),
                  "is reads considered single stranded. (Double strand by default)")
            ("single_read_length,l", po::value<int>(),
                                              "single read length")
            ("pair_read_length,i", po::value<std::vector<int> >()->multitoken(),
                                              "pair read length")
            ("jf_path,n", po::value<std::string>(),
                      "specify jf path from jellyfish")
            ("kmer_path,m", po::value<std::string>(),
                        "specify kmer path from jellyfish")
            ("SE_read_path,s", po::value<std::string>(),
                        "specify single ended read path ")
            ("PE_read_path,p", po::value<std::vector<std::string> >()->multitoken(),
                                              "specify pair read path")
            ("reference,r", po::value<std::string>(),
                                        "reference path for the output to compare with")
            ("output_dir,o", po::value<std::string>(),
                                    "output dir path (required)")
        ;

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << "Usage: options_description [options]\n";
            std::cout << desc;
            exit(0);
        }

        if(vm.count("reference"))
        {
            lf.reference_seq_path =  vm["reference"].as<std::string>();
            std::cout << "reference is " << lf.reference_seq_path << std::endl;
        }

        if (vm.count("config"))
        {
            std::string setting_file_name = vm["config"].as< std::string >();
            boost::filesystem::path base_path = boost::filesystem::current_path();
            std::string base_path_str(base_path.c_str());
            //("/shannon_C_setting.txt");
            std::string setting_path = base_path_str + "/" + setting_file_name;

            if(!boost::filesystem::exists(setting_path))
            {
                std::cout << "\033[0;33m";
                std::cout << "setting file: No such file: " << setting_path << std::endl;
                std::cout << "\033[0m";
                exit(1);
            }
            parser_setting_file(setting_path, setting);
        }
        else
        {
            std::cout << "\033[0;33m";
            std::cout << "ERROR: setting file is required" << std::endl;
            std::cout << "\033[0m";
            exit(0);
        }

        setting.is_double_stranded = true;
        if (vm.count("kmer_length"))
            setting.kmer_length = vm["kmer_length"].as<int>();
        if(vm.count("single_strand"))
        {
            std::cout << "activate single strand" << std::endl;
            setting.is_double_stranded = false;
        }


        if (vm.count("single_read_length"))
        {
            setting.has_single = true;
            setting.single_read_length = vm["single_read_length"].as<int>();

            se_len = true;
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
            pe_len = true;
        }

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
            if(!is_abs_path(lf.input_kmer_path))
                convert_relative_path_to_abs(lf.input_kmer_path);
        }

        if (vm.count("SE_read_path"))
        {
            lf.input_read_path = vm["SE_read_path"].as< std::string >();
            if(!is_abs_path(lf.input_read_path))
                convert_relative_path_to_abs(lf.input_read_path);
            lf.has_single = true;
            se_path = true;
        }

        if (vm.count("PE_read_path")) {
            lf.has_pair = true;
            lf.input_read_path_1 = vm["PE_read_path"].as< std::vector<std::string> >().at(0);
            if(!is_abs_path(lf.input_read_path_1))
                convert_relative_path_to_abs(lf.input_read_path_1);

            lf.input_read_path_2 = vm["PE_read_path"].as< std::vector<std::string> >().at(1);
            if(!is_abs_path(lf.input_read_path_2))
                convert_relative_path_to_abs(lf.input_read_path_2);
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

        if (vm.count("output_dir"))
        {
            lf.output_path = vm["output_dir"].as<std::string >();
            if(!is_abs_path(lf.output_path))
                convert_relative_path_to_abs(lf.output_path);
            std::cout << "lf.output_path " << lf.output_path<< std::endl;
            add_directory_if_not_exist(lf.output_path);
            lf.reset_paths();
        }
        else
        {
            std::cout << "\033[0;33m";
            std::cout << "ERROR: output dir is required" << std::endl;
            std::cout << "\033[0m";
            exit(0);
        }
        return;
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        exit(1);
    }
}
