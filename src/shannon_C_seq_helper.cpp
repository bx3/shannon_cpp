#include "shannon_C_seq_helper.h"

void produce_summary_file(struct Local_files & lf)
{
    lf.summary_file_path = lf.output_path + "/summary.log";
    std::ofstream writer(lf.summary_file_path);

    writer << "start at  " <<  lf.start_time << std::endl;
    writer << "finish at " <<  lf.end_time << std::endl;
    writer << "Notes: " << std::endl;
    //writer << setting->notes << std::endl;
    writer << std::endl;

    uint64_t num_output_seq = get_num_seq(lf.reconstructed_seq_path);
    size_t output_seq_size = std::ceil(get_filesize(lf.reconstructed_seq_path)/1024.0/1024.0);
    uint64_t num_reference_seq = get_num_seq(lf.reference_seq_path);
    size_t reference_seq_size = std::ceil(get_filesize(lf.reference_seq_path)/1024.0/1024.0);
    std::vector<std::string> py_summaries;
    std::string num_correct_seq = get_py_eval_summary(lf, py_summaries);

    int num_single_seq, num_contig_seq, num_sf_seq;
    count_seq_types_num(lf.reconstructed_seq_path,
                    num_single_seq, num_contig_seq, num_sf_seq);
    int align_num_single_seq, align_num_contig_seq, align_num_sf_seq;
    int num_single_contribute, num_contig_contribute, num_sf_contribute;
    int total_align = eval_align_counts(lf.eval_path + "_align" ,
                align_num_single_seq, align_num_contig_seq, align_num_sf_seq,
                num_single_contribute, num_contig_contribute,
                num_sf_contribute);
    int total_contribute = num_single_contribute+
                            num_contig_contribute+num_sf_contribute;



    int one_l = 1;
    writer << "               " << "aligned" << "\t|\t" << "output" << "\t|\t" << "percentage right" <<std::endl;
    writer << "total seq      " << total_align << "\t|\t" << num_output_seq
           << " (" << total_contribute << ")" << "\t|\t"
           << total_contribute*100/(std::max((uint64_t)1, num_output_seq))<<std::endl;
    writer << "single seq     " << align_num_single_seq << "\t|\t" << num_single_seq
           << " (" << num_single_contribute << ")" << "\t|\t"
           << num_single_contribute*100/(std::max(one_l, num_single_seq)) <<std::endl;
    writer << "contig seq     " << align_num_contig_seq << "\t|\t" << num_contig_seq
           << " (" << num_contig_contribute << ")"
           <<  "\t|\t" << num_contig_contribute*100/(std::max(one_l, num_contig_seq))<<std::endl;
    writer << "sf seq         " << align_num_sf_seq << "\t|\t" << num_sf_seq
           << " (" << num_sf_contribute << ")"
           <<  "\t|\t" << num_sf_contribute*100/(std::max(one_l, num_sf_seq))<<std::endl;
    writer << "output_seq_size    " << output_seq_size << " MiB" << std::endl;
    writer << std::endl;

    writer << "num_reference_seq  " << num_reference_seq << std::endl;
    writer << "reference_seq_size " << reference_seq_size << " MiB" << std::endl;
    writer << std::endl;

    for(int i=0; i<py_summaries.size();i++ )
        writer << py_summaries[i] << std::endl;

    //std::string setting_str = get_setting_string(*setting);
    //writer << std::endl;
    //writer << setting_str << std::endl;
    writer.close();

    std::string print_summary_cmd("cat " + lf.summary_file_path);
    run_command(print_summary_cmd, false);
}


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

        std::string local_cp_eval_file = lf.eval_dir_path + "/reference.fasta";
        std::string local_cp_eval_file_cmd =
            "cp " + lf.reference_seq_path + " " + local_cp_eval_file;
        run_command(local_cp_eval_file_cmd.c_str(), true);

        std::string cmd = "python " + curr_path + "/analysis/compare_trans.py " +
                          local_cp_eval_file + " " +
                          lf.reconstructed_seq_path + " " +
                          lf.eval_path;
        std::cout << cmd << std::endl;
        system(cmd.c_str());

        std::string rm_local_eval_file_cmd = "rm " + local_cp_eval_file;
        run_command(rm_local_eval_file_cmd.c_str(), true);
        std::string rm_local_cp_eval_nospace_file_cmd = "rm " + local_cp_eval_file + "_nospace";
        run_command(rm_local_cp_eval_nospace_file_cmd.c_str(), true);
        produce_summary_file(lf);
    }
    else
    {
        std::cout << "reference not provided, no evaluation" << std::endl;
    }
}


void eval_reconstructed_seq(std::string reference_seq_path,
                            std::string reconstructed_seq_path,
                            std::string eval_dir_path,
                            std::string eval_path)
{
    if(system(NULL))
        ;
    else
        exit(EXIT_FAILURE);

    if(reference_seq_path.empty())
    {
        std::cout << "no reference available" << std::endl;
        return;
    }

    if(!reference_seq_path.empty() && exist_path(reference_seq_path))
    {
        std::string curr_path = "./";
        convert_relative_path_to_abs(curr_path);

        std::string local_cp_eval_file = eval_dir_path + "/reference.fasta";
        std::string local_cp_eval_file_cmd =
            "cp " + reference_seq_path + " " + local_cp_eval_file;
        run_command(local_cp_eval_file_cmd.c_str(), true);

        std::string cmd = "python " + curr_path + "/analysis/compare_trans.py " +
                          local_cp_eval_file + " " +
                          reconstructed_seq_path + " " +
                          eval_path;
        std::cout << cmd << std::endl;
        system(cmd.c_str());

        std::string rm_local_eval_file_cmd = "rm " + local_cp_eval_file;
        run_command(rm_local_eval_file_cmd.c_str(), true);
        std::string rm_local_cp_eval_nospace_file_cmd = "rm " + local_cp_eval_file + "_nospace";
        run_command(rm_local_cp_eval_nospace_file_cmd.c_str(), true);
    }
    else
    {
        std::cout << "reference not provided, no evaluation" << std::endl;
    }
}

void count_seq_types_num(std::string filename, int & num_single_seq,
                        int & num_contig_seq, int & num_sf_seq)
{
    num_single_seq = 0;
    num_contig_seq = 0;
    num_sf_seq = 0;

    std::ifstream reader(filename);
    std::string line;
    while(std::getline(reader, line))
    {
        if(line[0] == '>')
        {
            char start_letter = line[1];
            if(start_letter=='s')
                num_single_seq++;
            else if(start_letter == 'C' || start_letter == 'S') // S for contigs from python Shannon
                num_contig_seq++;
            else if(start_letter =='c')
                num_sf_seq ++;
            else
            {
                std::cerr << "unknown type seq" << std::endl;
                exit(0);
            }
        }
    }
}

std::string get_py_eval_summary(Local_files & lf, std::vector<std::string> & py_summaries)
{
    std::ifstream reader(lf.eval_path);
    std::string line;
    std::string last_line;
    while(std::getline(reader, line))
    {
        if(line[0] == '#')
        {
            py_summaries.push_back(line);
        }
    }
    reader.close();
    return "";//line.substr(i);
}

int eval_align_counts(std::string filename,
        int & num_single_seq, int & num_contig_seq, int & num_sf_seq,
        int & num_single_contribute, int & num_contig_contribute,
        int & num_sf_contribute)
{
    num_single_seq = 0;
    num_contig_seq = 0;
    num_sf_seq = 0;

    std::ifstream reader(filename);
    std::string ref_name, recon_name, rest;
    std::set<std::string> single_contributing_seq;
    std::set<std::string> comp_contributing_seq;
    std::set<std::string> contig_contributing_seq;
    int total_align = 0;
    while(std::getline(reader, ref_name, '\t') &&
          std::getline(reader, recon_name, '\t') &&
          std::getline(reader, rest))
    {
        total_align++;
        char start_letter = recon_name[0];
        //std::cout << "recon_name" << recon_name << std::endl;
        //std::cout << "start_letter" << start_letter << std::endl;
        if(start_letter=='s')
        {
            num_single_seq++;
            single_contributing_seq.insert(recon_name);
        }
        else if(start_letter == 'C')
        {
            num_contig_seq++;
            contig_contributing_seq.insert(recon_name);
        }
        else if(start_letter =='c')
        {
            num_sf_seq ++;
            comp_contributing_seq.insert(recon_name);
        }
        else
        {
            shc_log_error("unknown type seq");
            exit(0);
        }
    }
    num_single_contribute = single_contributing_seq.size();
    num_contig_contribute = contig_contributing_seq.size();
    num_sf_contribute = comp_contributing_seq.size();
    return total_align;
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


void parse_jf_info(Local_files & lf, JF_stats & jf_stats)
{
    std::string cmd("jellyfish stats ");
    std::string jf_info_file(lf.output_path+ "/jf_stats");

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

    if(jf_stats.total_num == 0)
    {
        std::cerr << "[Error] jf stats file error, check if jellyfish stats run properly, "
                  << "total number of kmer shown to be 0 " << std::endl;
        exit(0);
    }


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

            std::cout << " Finish jellyfish" << std::endl;
            shc_log_info(shc_logname, "finish jellyfish with kmer length %d\n",
                                                setting.kmer_length);
        }
        else
        {
            std::string msg("Use existing jellyfish dict at path\n\t");
            msg += lf.input_kmer_path + "\n";
            msg += "\tdelete this file if input reads have changed.\n\n";
            msg += "\tUse jellyfish jf file at path\n\t";
            msg += lf.input_jf_path + "\n";
            msg += "\tdelete this file if input reads have changed.\n\n";
            print_important_notice(msg);
        }
    }
}

int get_kmer_length_from_kmer_file(std::string kmer_path)
{
    std::ifstream file_reader(kmer_path.c_str());
    std::string count_str, kmer_base;
    int kmer_length = -1;

    while(  std::getline(file_reader, kmer_base, '\t') &&
            std::getline(file_reader, count_str)   )
    {
        kmer_length = kmer_base.size();
        break;
    }
    assert(kmer_length > 0 && kmer_length <=32);
    return kmer_length;

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
