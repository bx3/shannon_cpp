#include "shannon_C_seq_helper.h"

void eval_reconstructed_seq(Local_files & lf)
{
    if(system(NULL))
        ;
    else
        exit(EXIT_FAILURE);

    if(lf.reference_seq_path.empty())
    {
        std::cout << "no reference available" << std::endl;
        return;
    }

    if(!lf.reference_seq_path.empty() && exist_path(lf.reference_seq_path))
    {
        std::string curr_path = "./";
        convert_relative_path_to_abs(curr_path);

        std::string cmd = "python " + curr_path + "/analysis/compare_trans.py " +
                          lf.reference_seq_path + " " +
                          lf.reconstructed_seq_path + " " +
                          lf.eval_path;
        std::cout << cmd << std::endl;
        system(cmd.c_str());
    }
    else
    {
        std::cout << "reference not provided, no evaluation" << std::endl;
    }
}

void check_read_type_consistency(Shannon_C_setting & setting)
{
    Local_files & lf = setting.local_files;
    if(setting.has_single)
    {
        assert(setting.single_read_length > 0);
        assert(!lf.input_read_path.empty());
    }
    if(setting.has_pair)
    {
        assert(setting.pair_1_read_length > 0);
        assert(setting.pair_2_read_length > 0);
        assert(!lf.input_read_path_1.empty());
        assert(!lf.input_read_path_2.empty());
    }
}

void
modify_setting_with_command(Local_files & lf, boost::program_options::variables_map & vm)
{
    if (vm.count("kmer_path"))
    {
        lf.input_kmer_path = vm["kmer_path"].as<std::string>();
        if(!is_abs_path(lf.input_kmer_path))
            convert_relative_path_to_abs(lf.input_kmer_path);
    }
    else
    {
        lf.input_kmer_path = "";
    }

    if (vm.count("SE_read_path"))
    {
        lf.input_read_path = vm["SE_read_path"].as< std::string >();
        if(!is_abs_path(lf.input_read_path))
            convert_relative_path_to_abs(lf.input_read_path);
        lf.has_single = true;
        std::cout << "SE read path are: "
             << lf.input_read_path << "\n";
    }

    if (vm.count("PE_read_path")) {
        lf.has_pair = true;
        lf.input_read_path_1 = vm["PE_read_path"].as< std::vector<std::string> >().at(0);
        if(!is_abs_path(lf.input_read_path_1))
            convert_relative_path_to_abs(lf.input_read_path_1);

        lf.input_read_path_2 = vm["PE_read_path"].as< std::vector<std::string> >().at(1);
        if(!is_abs_path(lf.input_read_path_2))
            convert_relative_path_to_abs(lf.input_read_path_2);
    }
    if (vm.count("jf_path")) {
        lf.input_jf_path = vm["jf_path"].as<std::string >();
        if(!is_abs_path(lf.input_jf_path))
            convert_relative_path_to_abs(lf.input_jf_path);
    }
    else
    {
        lf.input_jf_path = "";
    }

    if (vm.count("output_path")) {
        lf.output_path = vm["output_path"].as<std::string >();
        if(!is_abs_path(lf.output_path))
            convert_relative_path_to_abs(lf.output_path);
    }
    if (vm.count("reference")) {
        lf.reference_seq_path = vm["reference"].as<std::string >();
        if(!is_abs_path(lf.reference_seq_path))
            convert_relative_path_to_abs(lf.reference_seq_path);
    }
    lf.reset_paths();
}

int parse_main_command_line(int ac, char** av, Shannon_C_setting & setting)
{
    namespace po = boost::program_options;
    try {
        std::string config_path;
        int kmer_length;
        bool has_single = false;;
        bool has_pair = false;
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("config,j", po::value<std::string>(),
                  "json config file as default")
            ("kmer_length,k", po::value<int>()->default_value(25),
                  "specify kmer length.")
            ("kmer_path,m", po::value<std::string>(),
                        "specify kmer path from jellyfish")
            ("jf_path,f", po::value<std::string>(),
                                    "specify jellyfish jf path")
            ("SE_read_path,s", po::value<std::string>(),
                        "specify single ended read path ")
            ("PE_read_path,p", po::value< std::vector<std::string> >()->multitoken(),
                                    "specify paired ended read path ")
            ("output_path,o", po::value<std::string>(),
                                    "output dir path")
            ("reference,r", po::value<std::string>(),
                                            "reference dir path")
            ("double_stranded,d", "reference dir path")
            ("single_read_length,l", po::value<int>(),
                                            "single read length")
            ("pair_read_length,i", po::value<std::vector<int> >()->multitoken(),
                                            "pair read length")
            ("use_set,u", "use set to filter kmer")
            ("num_parallel,t", po::value<int>(), "number of threads for multi-graphs")
            ("load_factor,z", po::value<double>(), "load factor for kmer map")
        ;

        po::positional_options_description p;
        p.add("s", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac, av).
                  options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << "Usage: options_description [options]\n";
            std::cout << desc;
            return 1;
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
                shc_log_error("setting file: No such file: %s\n",
                                                setting_path.c_str());
                exit(1);
            }
            parser_setting_file(setting_path, setting);
        }
        else
        {
            std::cerr << "ERROR: a config file is required, use -j" << std::endl;
            return 1;
        }

        if (vm.count("double_stranded"))
        {
            setting.is_double_stranded = true;
        }
        else
        {
            setting.is_double_stranded = false;
        }

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

        if(vm.count("use_set"))
        {
            setting.dup_setting.is_use_set = true;
        }
        if(vm.count("num_parallel"))
        {
            setting.multi_graph_setup.num_parallel =
                          vm["num_parallel"].as<int>();
        }
        if(vm.count("load_factor"))
        {
            setting.dup_setting.load_factor = vm["load_factor"].as<double>();
        }

        modify_setting_with_command(setting.local_files, vm);

        return 0;
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        return 1;
    }
}

void parse_jf_info(Local_files & lf, JF_stats & jf_stats)
{
    std::string cmd("jellyfish stats ");
    std::string jf_info_file(lf.output_path+ "/jf_stats");

    std::cout << "lf.input_jf_path " << lf.input_jf_path << std::endl;
    std::cout << "jf_info_file     " << jf_info_file << std::endl;

    if(!exist_path(jf_info_file) || is_file_empty(jf_info_file))
    {
        cmd += (lf.input_jf_path + " > " + jf_info_file);
        run_command(cmd, true);
    }

    std::ifstream reader(jf_info_file);
    std::string name, num;
    while(std::getline(reader, name, ':') && std::getline(reader, num))
    {
        boost::trim(num);
        //std::cout << "num " << num << std::endl;

        if(name == "Unique")
            jf_stats.uniq_num = boost::lexical_cast<uint64_t>(num);
        else if (name == "Distinct")
            jf_stats.dist_num = boost::lexical_cast<uint64_t>(num);
        else if (name == "Total")
            jf_stats.total_num = boost::lexical_cast<uint64_t>(num);
        else if (name == "Max_count")
            jf_stats.max_count = boost::lexical_cast<uint64_t>(num);
        else
            std::cout << "read something strange" << name << std::endl;
    }

    std::cout << "uniq_num  " << jf_stats.uniq_num << "\n"
              << "dist_num  " << jf_stats.dist_num << "\n"
              << "total_num " << jf_stats.total_num << "\n"
              << "max_count " << jf_stats.max_count << std::endl;
}

void run_command(std::string cmd, bool print_cmd)
{
    if(system(NULL))
        ;
    else
        exit(EXIT_FAILURE);
    if(print_cmd)
    {
        std::cout << "\033[0;33m";
        std::cout << "run command" << std::endl;
        std::cout << cmd << std::endl;
        std::cout << "\033[0m" << std::endl << std::endl;
    }
    system(cmd.c_str());
}

void fasta_file_validator(std::string path)
{
    std::ifstream reader(path);
    std::string line;
    std::string next_line;
    int offset = 0;
    bool is_set_offset = false;
    int i = 1;
    int j = 0 ;

    while(std::getline(reader, line))
    {
        if(line.empty())
            break;
        i++;
        if(!is_set_offset && line.at(0) != '>')
        {
           is_set_offset = true;
           offset = i-1;
           std::cout << "offset = " << offset << std::endl;
        }
        else
        {
            j = i - offset;

            if(j%2 == 0 ) //
            {
                std::getline(reader, next_line);
                if(next_line.empty())
                    break;
                i++;
                char letter = next_line.at(0);
                if(letter == 'A' || letter == 'T' || letter == 'C' || letter == 'G')
                {
                    ;
                }
                else
                {
                    shc_log_error("seq at line %d has special %c \n", i, letter);
                    exit(1);
                }
            }
            else
            {
                shc_log_error("i=%d, j=%d : %s\n", i, j, line.c_str());
                exit(1);
            }
        }
    }
    std::cout << "it is valid, has read " << i << " line " << std::endl;
}

void run_jellyfish(Shannon_C_setting & setting)
{
    Local_files & lf = setting.local_files;

    if (lf.input_kmer_path.empty() )
    {
        lf.input_kmer_path = lf.algo_input + "/kmer.dict";
        lf.input_jf_path = lf.algo_input + "/kmer.jf";

        if(!exist_path(lf.input_kmer_path))
        {
            std::cout << "Run jellyfish with kmer length" << ((uint16_t)setting.kmer_length)
                      << std::endl;
            shc_log_info(shc_logname, "run jellyfish with kmer length %d\n",
                                                setting.kmer_length);

            std::string cmd_count = "jellyfish count -t 32 -o " + lf.input_jf_path +
                              " -m " + std::to_string(setting.kmer_length) +
                              " -s 200000000 ";

            std::string input_reads;
            if(setting.has_single)
            {
                input_reads += " " + lf.input_read_path;
            }

            if(setting.has_pair)
            {
                input_reads += " " + lf.input_read_path_1 + " " + lf.input_read_path_2;
            }

            if(setting.is_double_stranded)
                cmd_count += (std::string(" -C ") + input_reads);
            else
                cmd_count += input_reads;



            std::string cmd_dump = "jellyfish dump -ct -o " + lf.input_kmer_path +
                                " " + lf.input_jf_path;



            run_command(cmd_count, true);
            run_command(cmd_dump, true);

            std::cout << " Finish jellyfish: " << std::endl;
            shc_log_info(shc_logname, "finish jellyfish with kmer length %d\n",
                                                setting.kmer_length);
        }
    }
}

/*
int lock_memory(char   *addr,
            size_t  size)
{
  unsigned long    page_offset, page_size;

  page_size = sysconf(_SC_PAGE_SIZE);
  page_offset = (unsigned long) addr % page_size;

  addr -= page_offset;  //Adjust addr to page boundary
  size += page_offset;  //Adjust size with page_offset

  return ( mlock(addr, size) );  // Lock the memory
}

int unlock_memory(char   *addr,
              size_t  size)
{
  unsigned long    page_offset, page_size;

  page_size = sysconf(_SC_PAGE_SIZE);
  page_offset = (unsigned long) addr % page_size;

  addr -= page_offset;  // Adjust addr to page boundary
  size += page_offset;  // Adjust size with page_offset

  return ( munlock(addr, size) );  /* Unlock the memory
}
*/
