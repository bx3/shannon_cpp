#include "shannon_C_seq_helper.h"
#include <boost/tokenizer.hpp>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

uint64_t get_mem(pid_t id)
{
    FILE *in;
    char buff[512];

    char mypid[6];
    sprintf(mypid, "%d", id);
    std::string pid_string(mypid);

    std::string command = "ps -p " + pid_string + " -o rss | tail -n 2";
    if(!(in = popen(command.c_str(), "r"))){
        _exit(1);
    }
    while(fgets(buff, sizeof(buff), in)!=NULL){
        //cout << buff;
    }
    pclose(in);

    std::string str(buff);
    str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    uint64_t RSS = boost::lexical_cast<uint64_t>(str) ;
    return RSS;
}

uint64_t getMemoryUsage(pid_t id)
{
    std::string mem_log_file = "/proc/" + std::to_string(id)+ "/status";
    if(exist_path(mem_log_file))
    {
        FILE* file = fopen(mem_log_file.c_str(), "r");
        uint64_t result = 0;
        char line[128];

        while (fgets(line, 128, file) != NULL){
            if (strncmp(line, "VmRSS:", 6) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result*1024;  //return byte
    }
    else
        return 0;
}


void get_mem_statistics(int num_parallel, Mem_profiler & mp)
{
    std::string get_parallel_max_cmd = mp.parallel_dir;
    std::vector<std::string> lines(num_parallel);

    std::vector<std::shared_ptr<std::ifstream> > process_log_readers;
    for (int i=0; i<num_parallel; i++)
    {
        std::string process_log_path = mp.parallel_path_prefix + std::to_string(i) + ".ps.log";
        std::shared_ptr<std::ifstream> process_reader(new std::ifstream);

        process_reader->open(process_log_path.c_str());

        std::getline(*process_reader, lines[i]);

        process_log_readers.push_back(process_reader);

    }

    uint64_t max_RSS_sum = 0;
    uint64_t line_num = 1;
    uint64_t max_line = 0;

    typedef boost::tokenizer<boost::char_separator<char> >  tokenizer;
    boost::char_separator<char> sep(" ");

    while(true)
    {

        uint64_t RSS_sum = 0;
        for(int i=0; i<num_parallel; i++)
        {
            std::ifstream & reader = *(process_log_readers[i]);
            if(!std::getline(reader, lines[i]).eof())
            {

                tokenizer tok(lines[i], sep);
                int j=0;
                for(tokenizer::iterator it=tok.begin(); it!=tok.end();++it)
                {
                    if (j++ == 6)
                    {
                        uint64_t RSS = boost::lexical_cast<uint64_t>(*it);
                        //std::cout << "RSS " << RSS << std::endl;
                        RSS_sum +=RSS;
                    }
                }
            }
        }
        //std::cout << "RSS_sum "<< RSS_sum << std::endl;
        //std::cout << "line_num "<< line_num << std::endl;

        if(max_RSS_sum < RSS_sum)
        {
            max_RSS_sum = RSS_sum;
            max_line = line_num;
        }
        line_num++;
        if(RSS_sum == 0)
            break;
    }

    for (int i=0; i<num_parallel; i++)
    {
        (process_log_readers[i])->close();
    }


    std::ofstream writer(mp.mem_summary);
    writer << "paralle max_RSS_sum "<< max_RSS_sum/1024.0 << " MB <=> " << max_RSS_sum/1024.0/1014.0 << " GB" <<  std::endl;
    writer << "paralle max line    "<< max_line << std::endl;


    uint64_t main_RSS = get_max_RSS(mp.main_log_path+ ".ps.log");
    uint64_t jelly_RSS = get_max_RSS(mp.jelly_log_path+ ".ps.log");
    uint64_t sort_RSS = get_max_RSS(mp.sort_log_path+ ".ps.log");

    writer << "main_RSS  " << main_RSS/1024.0 << " MB <=> " << main_RSS/1024.0/1014.0 << " GB" << std::endl;
    writer << "jelly_RSS " << jelly_RSS/1024.0 << " MB <=> " << jelly_RSS/1024.0/1014.0 << " GB"<< std::endl;
    writer << "sort_RSS  " << sort_RSS/1024.0 << " MB <=> " << sort_RSS/1024.0/1014.0 << " GB" << std::endl;
    writer.close();
}



uint64_t get_max_RSS(std::string filename)
{
    std::ifstream reader(filename);
    std::string line;
    std::getline(reader, line);
    typedef boost::tokenizer<boost::char_separator<char> >  tokenizer;
    boost::char_separator<char> sep(" ");
    uint64_t RSS_max = 0;
    while(std::getline(reader, line))
    {

        tokenizer tok(line, sep);
        int j=0;
        for(tokenizer::iterator it=tok.begin(); it!=tok.end();++it)
        {
            if (j++ == 6)
            {
                uint64_t RSS = boost::lexical_cast<uint64_t>(*it);
                //std::cout << "RSS " << RSS << std::endl;
                if(RSS_max < RSS)
                    RSS_max = RSS;
            }
        }
    }
    reader.close();
    return RSS_max;
}

pid_t fork_mem_profiler(pid_t running_pid, std::string output_dir)
{
    pid_t profiler_pid;
    if( (profiler_pid = fork()) < 0)
    {
            printf("fork error");
            _exit(1);
    }
    else if (profiler_pid == 0)
    {
        std::string pid_str = std::to_string(running_pid);
        std::string command = "./syrupy.py";
        if (execlp(command.c_str(), command.c_str(), "-q", "-p",
            pid_str.c_str(), "-i", "10", output_dir.c_str(), (char *)0) < 0)
        {
                printf("execlp error");
                _exit(1);
        }
    }
    return profiler_pid;
}

void join_mem_profiler(pid_t profiler_pid)
{
    if (waitpid(profiler_pid, NULL, 0) < 0)
    {
            printf("wait profiler error");
            _exit(1);
    }
}

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
                num_sf_contribute, lf);
    int total_contribute = num_single_contribute+
                            num_contig_contribute+num_sf_contribute;

    writer << "number inside () denotes the number of unique aligned shannon output" << std::endl;
    writer << "so multiple ref transcripts might align to a single reconstructed transcript" << std::endl;

    int one_l = 1;
    writer << "               " << "aligned" << "\t|\t" << "sh output" << "\t|\t" << "percentage correct" <<std::endl;
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
        int & num_sf_contribute, struct Local_files & lf)
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
    std::ofstream comp_name_writer(lf.eval_dir_path + "/comp_align.fasta");
    for(std::set<std::string>::iterator it=comp_contributing_seq.begin();
                                    it!=comp_contributing_seq.end(); it++)
    {
        comp_name_writer << (*it) << std::endl;
    }

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

void print_yellow_cmd(std::string cmd)
{
    std::cout << "\033[0;33m";
    std::cout << "run command" << std::endl;
    std::cout << cmd << std::endl;
    std::cout << "\033[0m" << std::endl << std::endl;
}

bool run_command(std::string cmd, bool print_cmd)
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
    if(system(cmd.c_str())!=0)
        return false;
    return true;
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

// Assume the validity of fasta only check either fasta or fastq
bool is_fasta_file(std::string filename)
{
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(filename.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        if(seq->qual.l)
        {
            kseq_destroy(seq);
            gzclose(fp);
            return false;
        }
        else
        {
            kseq_destroy(seq);
            gzclose(fp);
            return true;
        }
    }

    std::cout << "\033[1;31m";
    std::cout << "Input is neither a fasta or fastq file" << std::endl;
    std::cout << "\033[0m" << std::endl;
    _exit(1);
}


bool is_gz_file(std::string filename)
{
    int len = filename.size();
    if (filename[len - 1] == 'z' && filename[len - 2] == 'g')
        return true;
    else
        return false;
}

void transform_to_fasta(
    std::string * setting_read_path,
    std::string output_name,
    std::string & raw_input_read_path,
    Local_files & lf)
{
    gzFile fp;
	kseq_t *seq;
	int l;
    fp = gzopen(raw_input_read_path.c_str(), "r");
    *setting_read_path = lf.algo_input + "/" + output_name;
    std::ofstream outfile (*setting_read_path);

	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
        outfile << ">" << seq->name.s << std::endl;
        outfile << seq->seq.s << std::endl;
	}

	kseq_destroy(seq);
	gzclose(fp);
}

void rcorrector_GetFileName( char *in, char *out )
{
	int i, j ;
	int len = (int)strlen( in ) ;
	for ( i = len ; i >= 0 && in[i] != '.' && in[i] != '/' ; --i )
		;
	if ( i >= 0 && !strcmp( &in[i], ".gz" ) )
	{
		int tmp = i ;
		for ( i = i - 1 ; i >= 0 && in[i] != '.' && in[i] != '/' ; --i )
			;
		in[tmp] = '\0' ;
		if ( i >= 0 && ( !strcmp( &in[i], ".fastq" ) || !strcmp( &in[i], ".fasta" ) ||
			!strcmp( &in[i], ".fq" ) || !strcmp( &in[i], ".fa" ) ) )
		{
			;
		}
		else
		{
			i = tmp ;
		}
		in[tmp] = '.' ;
	}

	for ( j = len ; j >= 0 && in[j] != '/' ; --j )
		;
	if ( i >= 0 && in[i] == '.' )
	{
		in[i] = '\0' ;
		strcpy( out, in + j + 1 ) ;
		in[i] = '.' ;
	}
	else
	{
		strcpy( out, in + j + 1 ) ;
	}
}

std::string get_filename_rcorrector(std::string & input)
{
    char fileName[1024];
    char in[1024];
    strcpy(in, input.c_str());
    rcorrector_GetFileName(in, fileName);
    std::string filename(fileName);
    return filename;
}

void run_pre_error_correct(Shannon_C_setting & setting)
{
    std::cout << "run_pre_error_correct" << std::endl;

    Local_files & lf = setting.local_files;
    std::string num_thread_str = std::to_string(setting.num_parallel);

    int num_field = 15;
    char **argv = (char**) malloc(num_field*sizeof(char*));
    for(int i=0;i<num_field;i++){
        argv[i] = (char*) malloc(1000*sizeof(char));
    }

    std::string rcorrector_path = lf.shannon_env_path +  "/Rcorrector/run_rcorrector.pl";
    if (!exist_path(rcorrector_path))
    {
        rcorrector_path = std::string("run_rcorrector.pl");
    }

    std::cout << "rcorrector_path " << rcorrector_path << std::endl;
    int arg_n = 0;
    //strcpy(argv[arg_n++], "perl");
    strcpy(argv[arg_n++], rcorrector_path.c_str());
    strcpy(argv[arg_n++], "-t");
    strcpy(argv[arg_n++], num_thread_str.c_str());
    strcpy(argv[arg_n++], "-od");
    strcpy(argv[arg_n++], lf.algo_input.c_str());


    if(setting.has_single)
    {

        std::string filename = get_filename_rcorrector(lf.input_read_path);

        if(is_fasta_file(lf.input_read_path))
            lf.pre_corrected_read_path = lf.algo_input + "/" + filename + ".cor.fa";
        else
            lf.pre_corrected_read_path = lf.algo_input + "/" + filename + ".cor.fq";

        if(is_gz_file(lf.input_read_path))
            lf.pre_corrected_read_path += ".gz";

        strcpy(argv[arg_n++], "-r");
        strcpy(argv[arg_n++], lf.input_read_path.c_str());
    }

    if(setting.has_pair)
    {
        std::string filename1 = get_filename_rcorrector(lf.input_read_path_1);
        std::string filename2 = get_filename_rcorrector(lf.input_read_path_2);
        if (is_fasta_file(lf.input_read_path_1))
            lf.pre_corrected_read_path_1 = lf.algo_input + "/" + filename1 +".cor.fa";
        else
            lf.pre_corrected_read_path_1 = lf.algo_input + "/" + filename1 +".cor.fq";
        if(is_gz_file(lf.input_read_path_1))
            lf.pre_corrected_read_path_1 += ".gz";


        if (is_fasta_file(lf.input_read_path_2))
            lf.pre_corrected_read_path_2 = lf.algo_input + "/" + filename2 +".cor.fa";
        else
            lf.pre_corrected_read_path_2 = lf.algo_input + "/" + filename2 +".cor.fq";
        if(is_gz_file(lf.input_read_path_2))
            lf.pre_corrected_read_path_2 += ".gz";


        strcpy(argv[arg_n++], "-1");
        strcpy(argv[arg_n++], lf.input_read_path_1.c_str());
        strcpy(argv[arg_n++], "-2");
        strcpy(argv[arg_n++], lf.input_read_path_2.c_str());
    }

    if(setting.has_single && setting.has_pair)
    {
        if(exist_path(lf.pre_corrected_read_path) &&
           exist_path(lf.pre_corrected_read_path_1) &&
           exist_path(lf.pre_corrected_read_path_2)   )
           {
               lf.input_read_path = lf.pre_corrected_read_path;
               lf.input_read_path_1 = lf.pre_corrected_read_path_1;
               lf.input_read_path_2 = lf.pre_corrected_read_path_2;
               return;
           }
    }
    else if(setting.has_single)
    {
        if (exist_path(lf.pre_corrected_read_path))
        {
            lf.input_read_path = lf.pre_corrected_read_path;
            return;
        }
    }
    else
    {
        if (exist_path(lf.pre_corrected_read_path_1) &&
            exist_path(lf.pre_corrected_read_path_2)   )
        {
            lf.input_read_path_1 = lf.pre_corrected_read_path_1;
            lf.input_read_path_2 = lf.pre_corrected_read_path_2;
            return;
        }
    }


    argv[arg_n++] = NULL;

    std::string cmd_count;
    for(int i=0;i<arg_n-1;i++)
    {
        cmd_count += std::string(argv[i]) + " ";
    }

    //print_yellow_cmd(cmd_count);
    if(!run_command(cmd_count, true))
        exit(EXIT_FAILURE);
    //if(execvp("perl", argv))
    //{
    //    printf("execlp error");
    //    _exit(1);
    //}

    if(setting.has_single)
    {
        lf.input_read_path = lf.pre_corrected_read_path;
    }

    if(setting.has_pair)
    {
        lf.input_read_path_1 = lf.pre_corrected_read_path_1;
        lf.input_read_path_2 = lf.pre_corrected_read_path_2;
    }

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
            std::string num_thread_str = std::to_string(setting.num_parallel);
            std::string kmer_length_str = std::to_string(setting.kmer_length);

            int num_field = 15;
            char **argv = (char**) malloc(num_field*sizeof(char*));
            for(int i=0;i<num_field;i++){
                argv[i] = (char*) malloc(1000*sizeof(char));
            }

            int arg_n = 0;
            strcpy(argv[arg_n++], "jellyfish");
            strcpy(argv[arg_n++], "count");
            strcpy(argv[arg_n++], "-t");
            strcpy(argv[arg_n++], num_thread_str.c_str());
            strcpy(argv[arg_n++], "-o");
            strcpy(argv[arg_n++], lf.input_jf_path.c_str());
            strcpy(argv[arg_n++], "-m");
            strcpy(argv[arg_n++], kmer_length_str.c_str());
            strcpy(argv[arg_n++], "-s");
            strcpy(argv[arg_n++], "200000000");



            //std::string cmd_count = "jellyfish count -t 16 -o " + lf.input_jf_path +
            //                  " -m " + std::to_string(setting.kmer_length) +
            //                  " -s 200000000 ";
            if(setting.is_double_stranded)
            {
                //cmd_count += (std::string(" -C "));
                strcpy(argv[arg_n++], "-C");
            }

            std::string input_reads;
            if(setting.has_single)
            {
                //input_reads += lf.input_read_path;
                strcpy(argv[arg_n++], lf.input_read_path.c_str());
                argv[arg_n++] =  NULL;
            }

            if(setting.has_pair)
            {
                //input_reads += lf.input_read_path_1 + " " + lf.input_read_path_2;
                strcpy(argv[arg_n++], lf.input_read_path_1.c_str());
                strcpy(argv[arg_n++], lf.input_read_path_2.c_str());
                argv[arg_n++] = NULL;
            }


            // cmd_count += input_reads;

            std::string cmd_dump = "jellyfish dump -ct -o " + lf.input_kmer_path +
                                " " + lf.input_jf_path;

            std::string cmd_count;
            for(int i=0;i<arg_n-1;i++)
            {
                cmd_count += std::string(argv[i]) + " ";
            }

            pid_t pid;

            if ((pid = fork()) < 0)
            {
                printf("fork error");
                _exit(1);
            }
            else if (pid == 0)
            {
                print_yellow_cmd(cmd_count);
                if(execvp("jellyfish", argv))
                {
                    printf("execlp error");
                    _exit(1);
                }

                _exit(0);
            }

            pid_t profiler_pid;
            if( (profiler_pid = fork()) < 0)
            {
                    printf("fork error");
                    _exit(1);
            }
            else if (profiler_pid == 0)
            {
                    printf("jelly pid %d\n", pid);
                    std::string pid_str = std::to_string(pid);

                    if (execlp("./syrupy.py", "./syrupy.py", "-q", "-p", pid_str.c_str(), "-i", "10",
                                setting.local_files.mem_profiler.jelly_log_path.c_str(),(char *)0) < 0)
                    {
                            printf("execlp error");
                            _exit(1);
                    }
                    _exit(0);
            }

            if (waitpid(pid, NULL, 0) < 0)
            {
                    printf("wait error");
                    _exit(1);
            }
            if (waitpid(profiler_pid, NULL, 0) < 0)
            {
                    printf("wait profiler error");
                    _exit(1);
            }

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

/*
if(setting.is_double_stranded)
{


    if (execlp("jellyfish", "jellyfish", "count", "-t", num_thread_str.c_str(),
                "-o", lf.input_jf_path.c_str(),
                "-m", kmer_length_str.c_str(),
                "-s", "200000000", "-C",  input_reads.c_str(), (char *)0) < 0)
    {
            printf("execlp error");
            _exit(1);
    }
}
else
{
    if (execlp("jellyfish", "jellyfish", "count", "-t", num_thread_str.c_str(),
                "-o", lf.input_jf_path.c_str(),
                "-m", kmer_length_str.c_str(),
                "-s", "200000000", input_reads.c_str(), (char *)0) < 0)
    {
            printf("execlp error");
            _exit(1);
    }
}
*/

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




void print_and_log_all_setting(Shannon_C_setting & setting)
{
    print_and_log_general_setting(setting);
    print_and_log_partition_setting(setting);
    print_and_log_local_file_system(&setting.local_files);
    print_and_log_multi_graph_setting(setting);
    print_and_log_mb_setting(setting);
    print_and_log_sf_setting(setting);
    print_and_log_find_rep_setting(setting);
}

void print_and_log_find_rep_setting(Shannon_C_setting & setting)
{
    std::cout << "\033[1;34m";
    std::cout << "find-rep         ***********" << std::endl;
    std::cout << "rmer_length is           : " << setting.rmer_length << std::endl;
    std::cout << "\033[0m" << std::endl;
    shc_log_info(shc_logname, "rmer_length is             : %d\n", setting.rmer_length);
}

void print_and_log_shannon_cmd_setting(Shannon_C_setting & setting)
{
    print_and_log_all_setting(setting);
}

void print_and_log_partition_cmd_setting(Shannon_C_setting & setting)
{
    print_and_log_general_setting(setting);
    print_and_log_partition_setting(setting);
}

void print_and_log_multi_graph_both_cmd_setting(Shannon_C_setting & setting)
{
    print_and_log_general_setting(setting);
    print_and_log_multi_graph_setting(setting);
    print_and_log_mb_setting(setting);
    print_and_log_sf_setting(setting);
}

void print_and_log_output_setting(Shannon_C_setting & setting)
{
    Local_files & local_files = setting.local_files;
    std::cout << "\033[1;33m";
    std::cout << "output save to           : " << local_files.output_path << std::endl;
    std::cout << "output_seq_min_len is    : " << setting.output_seq_min_len << std::endl;
    std::cout << "\033[0m" << std::endl;
    shc_log_info(shc_logname, "output save to             : %s\n", local_files.output_kmer_path.c_str());
    shc_log_info(shc_logname, "output_seq_min_len         : %d\n", setting.output_seq_min_len);
}

void print_and_log_ref_setting(Shannon_C_setting & setting)
{
    std::cout << "\033[1;33m";
    if(!setting.local_files.reference_seq_path.empty())
        std::cout << "reference path           : " << setting.local_files.reference_seq_path << std::endl;
    else
        std::cout << "reference path           : No ref available" << std::endl;
    std::cout << "\033[0m" << std::endl;

    if(!setting.local_files.reference_seq_path.empty())
        shc_log_info(shc_logname, "reference path           : %s\n", setting.local_files.reference_seq_path.c_str());
    else
        shc_log_info(shc_logname, "reference path           : No ref available\n");
}

void print_and_log_partition_setting(Shannon_C_setting & setting)
{
    std::cout << "\033[1;31m";  //34 blue, 31 red, 35 purple
    std::cout << "Duplicate correction     *********** " << std::endl;
    std::cout << "load_factor              : " << setting.dup_setting.load_factor << std::endl;
    std::cout << "rmer_length              : " << (uint16_t)setting.dup_setting.rmer_length << std::endl;
    std::cout << "min_count                : " << setting.dup_setting.min_count << std::endl;
    std::cout << "min_len                  : " << setting.dup_setting.min_len << std::endl;
    std::cout << "threshold                : " << setting.dup_setting.threshold << std::endl;
    std::cout << "is_use_set               : " << ((setting.dup_setting.is_use_set)?("Yes"):("No")) << std::endl;
    std::cout << "num_sort_thread          : " << ((setting.dup_setting.num_sort_thread)) << std::endl;
    std::cout << "sort_tmp_dir           : " << ((setting.dup_setting.sort_tmp_dir)) << std::endl;
    std::cout << std::endl;

    std::cout << "\033[1;31m";  //34 blue, 31 red, 35 purple
    std::cout << "Metis setup              *********** " << std::endl;
    std::cout << "use_multiple_partition   : " << ((setting.metis_setup.is_multiple_partition)?("Yes"):("No")) << std::endl;
    std::cout << "partition_size           : " << setting.metis_setup.partition_size << std::endl;
    std::cout << "non_partition_size       : " << setting.metis_setup.non_partition_size << std::endl;
    std::cout << "penalty                  : " << setting.metis_setup.penalty << std::endl;
    std::cout << "overload                 : " << setting.metis_setup.overload << std::endl;
    std::cout << std::endl;

    std::cout << "\033[1;31m";  //34 blue, 31 red, 35 purple
    std::cout << "Contig graph setup       *********** " << std::endl;
    std::cout << "num_feature                 : " << setting.contig_graph_setup.num_feature << std::endl;
    std::cout << "is_assign_best           : " << ((setting.contig_graph_setup.is_assign_best)?("Yes"):("No")) << std::endl;
    std::cout << "read_sampler_k           : " << setting.contig_graph_setup.read_sampler_k << std::endl;
    std::cout << "\033[0m" << std::endl;
    std::cout << std::endl;

    shc_log_info(shc_logname, "Duplicate correction     ***********\n");
    shc_log_info(shc_logname, "load_factor              : %f\n", setting.dup_setting.load_factor);
    shc_log_info(shc_logname, "rmer_length              : %d\n", (uint16_t)setting.dup_setting.rmer_length);
    shc_log_info(shc_logname, "min_count                : %d\n", setting.dup_setting.min_count);
    shc_log_info(shc_logname, "min_len                  : %d\n", setting.dup_setting.min_len);
    shc_log_info(shc_logname, "threshold                : %f\n", setting.dup_setting.threshold);
    shc_log_info(shc_logname, "is_use_set               : %s\n", ((setting.dup_setting.is_use_set)?("Yes"):("No")));
    shc_log_info(shc_logname, "num_sort_thread          : %d\n", ((setting.dup_setting.num_sort_thread)));
    shc_log_info(shc_logname, "\n");

    shc_log_info(shc_logname, "Metis setup              *********** \n");
    shc_log_info(shc_logname, "use_multiple_partition   : %s\n", ((setting.metis_setup.is_multiple_partition)?("Yes"):("No")));
    shc_log_info(shc_logname, "partition_size           : %d\n", setting.metis_setup.partition_size);
    shc_log_info(shc_logname, "non_partition_size       : %d\n", setting.metis_setup.non_partition_size);
    shc_log_info(shc_logname, "penalty `                : %d\n", setting.metis_setup.penalty);
    shc_log_info(shc_logname, "overload                 : %f\n", setting.metis_setup.overload);
    shc_log_info(shc_logname, "\n");

    shc_log_info(shc_logname, "Contig graph setup       *********** \n");
    shc_log_info(shc_logname, "num_feature                 : %d\n", setting.contig_graph_setup.num_feature);
    shc_log_info(shc_logname, "is_assign_best           : %s\n", ((setting.contig_graph_setup.is_assign_best)?("Yes"):("No")));
    shc_log_info(shc_logname, "read_sampler_k           : %d\n", setting.contig_graph_setup.read_sampler_k);

    shc_log_info(shc_logname, "\n");
}

void print_and_log_input_path_setting(Shannon_C_setting & setting)
{
    Local_files & local_files = setting.local_files;
    std::cout << "\033[1;33m";
    if(!local_files.input_kmer_path.empty())
    {
        std::cout << "single kmer input from   : " << local_files.input_kmer_path << std::endl;
    }
    if(local_files.has_single)
    {
        std::cout << "single read input from   : " << local_files.input_read_path << std::endl;
    }

    if(local_files.has_pair)
    {
        std::cout << "read pair1 input from    : " << local_files.input_read_path_1 << std::endl;
        std::cout << "read pair2 input from    : " << local_files.input_read_path_2 << std::endl;
    }
    std::cout << "\033[0m" << std::endl;


    if(!local_files.input_kmer_path.empty())
    {
        shc_log_info(shc_logname, "single kmer input from   : %s", local_files.input_kmer_path.c_str());
    }
    if(local_files.has_single)
    {
        shc_log_info(shc_logname, "single read input from   : %s", local_files.input_read_path.c_str());
    }

    if(local_files.has_pair)
    {
        shc_log_info(shc_logname, "read pair1 input from    : %s\n", local_files.input_read_path_1.c_str());
        shc_log_info(shc_logname, "read pair2 input from    : %s\n", local_files.input_read_path_2.c_str());
    }
    if(!local_files.input_kmer_path.empty())
    {
        shc_log_info(shc_logname, "single kmer input from   : %s\n", local_files.input_kmer_path.c_str());
    }
}

void print_and_log_read_length_setting(Shannon_C_setting & setting)
{
    std::cout << "\033[1;33m";
    std::cout << "has_single               : " << ((setting.has_single)?("Yes"):("No")) << std::endl;
    if(setting.has_single)
        std::cout << "single read length       : " << (setting.single_read_length) << std::endl;
    std::cout << "has_pair                 : " << ((setting.has_pair)?("Yes"):("No")) << std::endl;
    if(setting.has_pair)
    {
        std::cout << "pair 1 length            : " << (setting.pair_1_read_length) << std::endl;
        std::cout << "pair 2 length            : " << (setting.pair_2_read_length) << std::endl;
    }
    std::cout << "\033[0m";

    shc_log_info(shc_logname, "has_single               : %s\n", ((setting.has_single)?("Yes"):("No")));
    if(setting.has_single)
        shc_log_info(shc_logname, "single length            : %d\n", (setting.single_read_length));
    shc_log_info(shc_logname, "has_pair                 : %s\n", ((setting.has_pair)?("Yes"):("No")));
    if(setting.has_pair)
    {
        shc_log_info(shc_logname, "pair 1 length            : %d\n", (setting.pair_1_read_length));
        shc_log_info(shc_logname, "pair 2 length            : %d\n", (setting.pair_2_read_length));
    }
}

void print_and_log_kmer_strand_setting(Shannon_C_setting & setting)
{
    std::cout << "\033[1;33m";
    std::cout << "kmer_length is           : " << static_cast<uint16_t>(setting.kmer_length) << std::endl;
    std::cout << "double stranded          : " << ((setting.is_double_stranded)?("Yes"):("No")) << std::endl;
    std::cout << "\033[0m";
    shc_log_info(shc_logname, "kmer_length is           : %d\n", static_cast<uint16_t>(setting.kmer_length));
    shc_log_info(shc_logname, "double stranded          : %s\n", ((setting.is_double_stranded)?("Yes"):("No")));
}

void print_and_log_general_setting(Shannon_C_setting & setting)
{
    ;
}

void print_and_log_mb_setting(Shannon_C_setting & setting)
{
    std::cout << "\033[1;34m";
    std::cout << "Multi-bridge     ***********" << std::endl;
    std::cout << "max_hop_path             : " << setting.seq_graph_setup.max_hop_path << std::endl;
    std::cout << "\033[0m" << std::endl;
    shc_log_info(shc_logname, "Multi-bridge     ***********\n");
    shc_log_info(shc_logname, "max_hop_path             : %d\n", setting.seq_graph_setup.max_hop_path);
}

void print_and_log_sf_setting(Shannon_C_setting & setting)
{
    std::cout << "\033[1;34m";
    std::cout << "Sparse flow     ***********" << std::endl;
    std::cout << "sparse flow multiple test: " << setting.sparse_flow_setup.multiple_test << std::endl;
    std::cout << "\033[0m" << std::endl;
    shc_log_info(shc_logname, "Sparse flow     ***********\n");
    shc_log_info(shc_logname, "sparse flow multiple test: %d\n", setting.sparse_flow_setup.multiple_test);
}

void print_and_log_multi_graph_setting(Shannon_C_setting & setting)
{
    std::cout << "\033[1;34m";
    std::cout << "Multi-graph     ***********" << std::endl;
    std::cout << "num_parallel             : " << setting.num_parallel << std::endl;
    std::cout << "avail mem (byte)         : " << setting.avail_mem << std::endl;
    std::cout << "\033[0m" << std::endl;
    shc_log_info(shc_logname, "Multi-graph     ***********\n");
    shc_log_info(shc_logname, "num_parallel             : %d\n", setting.num_parallel);
    shc_log_info(shc_logname, "avail mem (byte          : %d\n", setting.avail_mem);
}

//https://stackoverflow.com/questions/2513505/how-to-get-available-memory-c-g/26639774
#ifdef __unix__ 
int64_t get_machine_physical_limit_mem()
{
    std::cout << "detect linux/unix system" << std::endl;
    int64_t pages = sysconf(_SC_PHYS_PAGES);
    int64_t page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}
#endif

#ifdef __APPLE__
int64_t get_machine_physical_limit_mem()
{
    std::cout << "detect Mac OS" << std::endl;
    size_t memorySize;
    sysctlbyname("hw.memsize",nullptr,&memorySize,nullptr,0);

    return static_cast<int64_t>(memorySize);
}
#endif

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

  return ( munlock(addr, size) );  // Unlock the memory
}
*/
