#include "local_file_structure.h"

bool convert_to_abs_and_check_exist(std::string & path)
{
    if(!is_abs_path(path))
        convert_relative_path_to_abs(path);
    if(!exist_path(path))
    {
        std::cout << "\033[0;33m";
        std::cout << "path not exist"  << std::endl;
        std::cout << path  << std::endl;
        std::cout << "\033[0m";
        return false;
    }
    return true;
}

void convert_relative_path_to_abs(std::string rel_path, std::string & abs_path)
{
    char resolved_path[PATH_MAX];
    realpath(rel_path.c_str(), resolved_path);
    std::string abs_out_path(resolved_path);
    abs_path = abs_out_path;
}

void convert_relative_path_to_abs(std::string & a_path)
{
    char resolved_path[PATH_MAX];
    realpath(a_path.c_str(), resolved_path);
    std::string abs_out_path(resolved_path);
    a_path = abs_out_path;
}

bool is_abs_path(std::string & a_path)
{
    return a_path.at(0) == '/';
}

void overwirte_a_file(std::string file_path)
{
    std::ofstream file_writer(file_path);
    file_writer.close();
}

void copy_file(std::string file_path, std::string copy_file_path)
{
    if(system(NULL))
        std::cout << "copy " <<file_path<< " to " << copy_file_path<< std::endl;
    else
        exit(EXIT_FAILURE);

    std::string cmd("cp ");
    cmd = cmd + file_path + " " + copy_file_path;
    system(cmd.c_str());
}

int count_num_files(std::string path)
{
    int counter=0;
    for(boost::filesystem::directory_iterator it(path);
        it != boost::filesystem::directory_iterator(); it++)
    {
        counter++;
    }
    return counter;
}

bool exist_path(std::string a_path)
{
    boost::filesystem::path a_path_boost(a_path);
    return boost::filesystem::exists(a_path_boost);
}

void add_or_overwrite_directory(std::string & dir, std::string & output_dir)
{
    if(!is_abs_path(dir))
        convert_relative_path_to_abs(dir);
    if(!is_abs_path(output_dir))
        convert_relative_path_to_abs(output_dir);
    std::cout << "add or overwrite dir " << dir << std::endl;
    int output_path_length = output_dir.size();
    std::string dir_prefix = dir.substr(0, output_path_length);
    if( dir.size()<=output_path_length || dir_prefix != output_dir)
    {
        std::cout << "cannot overwrite, since dir is not a subdir of output dir" << std::endl;
        std::cout << "new dir : " << dir << std::endl;
        std::cout << "output dir: " << output_dir << std::endl;
        exit(0);
    }

    boost::filesystem::path dir_boost_path(dir);
    if(boost::filesystem::exists(dir))
        empty_directory(dir_boost_path, output_dir);
    else
        add_directory(dir_boost_path);
}

void add_directory_if_not_exist(std::string & dir)
{

    boost::filesystem::path dir_boost_path(dir);
    if(!boost::filesystem::exists(dir))
    {
        std::cout << "add dir " << dir << std::endl;
        add_directory(dir_boost_path);
    }
    else
    {
        std::cout << "dir already exists, use the existing one : "
                  << dir << std::endl << std::endl;
    }

}

bool is_file_empty(std::string a_path)
{
    std::ifstream reader(a_path);
    bool is_empty = reader.peek() == std::ifstream::traits_type::eof();
    reader.close();
    return is_empty;
}

void empty_directory(boost::filesystem::path & dir_path, std::string output_dir)
{
    int output_path_length = output_dir.size();
    std::string dir = dir_path.string();
    std::string dir_prefix = dir.substr(0, output_path_length);
    if( dir.size()<=output_path_length || dir_prefix != output_dir)
    {
        std::cout << "cannot clean, since dir is not a subdir of output dir" << std::endl;
        std::cout << "clean dir : " << dir << std::endl;
        std::cout << "output dir: " << output_dir << std::endl;
        exit(0);
    }
    std::cout << "*********clean dir" << std::endl;
    boost::filesystem::remove_all(dir_path);
    add_directory(dir_path);
}

void add_directory(boost::filesystem::path & dir_path)
{
    if(boost::filesystem::create_directory(dir_path))
    {
#ifdef SHOW_CREATE_DIR
        std::cout << "success to add directory: " << dir_path << std::endl;
#endif
    }
    else
    {
#ifdef SHOW_CREATE_DIR
        std::cout << "fail to add directory: " << dir_path << std::endl;
#endif
    }
}

void remove_directory(boost::filesystem::path & dir_path)
{
    char letter;
    bool flag = true;
    std::cout << "Are you sure to delete the following dir: (Y/N) " << std::endl;
    std::cout << dir_path << std::endl;
    std::cin >> letter;
    while(flag)
    {
        if(letter=='Y' || letter=='y')
        {
            boost::filesystem::remove_all(dir_path);
            std::cout << "successfully delete the folder" << std::endl;
            return;
        }
        else if(letter=='N' || letter=='n')
        {
            return;
        }
        else
        {
            std::cout << "please enter Y or N" << std::endl;
            std::cin >> letter;
        }
    }
}

void remove_file(const boost::filesystem::path & file_path)
{
    boost::filesystem::remove(file_path);
}

size_t get_filesize(const std::string & filename)
{
    if(boost::filesystem::exists(filename))
    {
        std::ifstream file_reader(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
        return file_reader.tellg();
    }
    return 0;
}

size_t estimate_num_read(const std::string & filename, size_t read_length)
{
    std::string test_filename(filename+"test_read");
    std::ofstream test_file(test_filename.c_str());
    for (int i=0; i<TEST_NUM_READ ; i++)
    {
        for(int j=0; j<10; j++)
            test_file << "H";
        test_file << std::endl;
        for(int j=0; j<read_length; j++)
            test_file << "A";
        test_file << std::endl;
    }
    test_file.close();
    size_t byte_per_line = get_filesize(test_filename) / TEST_NUM_READ;
    size_t estimated_lines = static_cast<size_t>
                        (static_cast<double>(get_filesize(filename)) /
                         static_cast<double>(byte_per_line));
    boost::filesystem::path test_file_path(test_filename);
    remove_file(test_file_path);
    return estimated_lines;
}

size_t estimate_num_kmer(const std::string & filename, uint8_t kmer_length)
{
    std::string test_filename(filename+"test_kmer");
    std::ofstream test_file(test_filename.c_str());
    for (int i=0; i<TEST_NUM_LINE ; i++)
    {
        for(int j=0; j<kmer_length; j++)
            test_file << "A";
        test_file << "\t" << 50  << std::endl;
    }
    test_file.close();
    size_t byte_per_line = get_filesize(test_filename) / TEST_NUM_LINE;
    size_t estimated_lines = static_cast<size_t>
                        (static_cast<double>(get_filesize(filename)) /
                         static_cast<double>(byte_per_line)*SIZE_MULTIPLIER);
    boost::filesystem::path test_file_path(test_filename);
    remove_file(test_file_path);
    return estimated_lines;
}


void print_and_log_local_file_system(Local_files *local_files)
{
    std::cout << "\033[1;32m";  //34 blue, 31 red, 35 purple
    std::cout << "Input files               *********** " << std::endl;


    std::cout << "Output                    *********** " << std::endl;
    std::cout << "output save to           : " << local_files->output_path << std::endl;

    std::cout << "\033[0m" << std::endl;

    shc_log_info(shc_logname, "current work dir is      : %s\n", local_files->base_path.c_str());
    shc_log_info(shc_logname, "\n");
    shc_log_info(shc_logname, "Input                     *********** \n");
    shc_log_info(shc_logname, "If has single read       : %s\n", ((local_files->has_single)?"Yes":"No"));
    shc_log_info(shc_logname, "current work dir is      : %s\n", local_files->base_path.c_str());
    shc_log_info(shc_logname, "reference path           : %s\n", local_files->reference_seq_path.c_str());

    shc_log_info(shc_logname, "single kmer input from   : %s\n", local_files->input_kmer_path.c_str());


    shc_log_info(shc_logname, "Output                    *********** \n");
    shc_log_info(shc_logname, "kmer output save to        : %s\n", local_files->output_kmer_path.c_str());
    shc_log_info(shc_logname, "contig output saved to     : %s\n", local_files->output_contig_path.c_str());
    shc_log_info(shc_logname, "METIS components under dir : %s\n", local_files->output_comp_path.c_str());
    shc_log_info(shc_logname, "MB node/edge/path under    : %s\n", local_files->output_seq_graph_path.c_str());
    shc_log_info(shc_logname, "SF reconstructed seq under : %s\n", local_files->reconstructed_seq_path.c_str());

    shc_log_info(shc_logname, "\n");
}


uint64_t get_num_seq(std::string in_file)
{
    FILE *cmd;
    char result[1024];
    std::string grep_cmd("grep -c '>' ");
    grep_cmd += in_file;
    cmd = popen(grep_cmd.c_str(), "r");
    if (cmd == NULL) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    while (fgets(result, sizeof(result), cmd)) {
        ;//printf("%s", result);
    }
    pclose(cmd);
    uint64_t num_seq = atoi(result);
    //std::cout << num_seq << std::endl;
    return num_seq;
}

uint64_t get_num_kmer(std::string in_file)
{
    FILE *cmd;
    char result[1024];
    std::string wc_cmd("wc -l  ");
    wc_cmd += in_file;
    cmd = popen(wc_cmd.c_str(), "r");
    if (cmd == NULL) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    while (fgets(result, sizeof(result), cmd)) {
        ;//printf("%s", result);
    }
    pclose(cmd);
    uint64_t num_seq = atoi(result);
    //std::cout << num_seq << std::endl;
    return num_seq;
}


void replace_space_with_underscore(std::string & s)
{
    for(int i=0; i<s.size(); i++)
    {
        char & letter = s[i];
        if(letter == ' ')
            letter = '_';
    }
}

std::string get_current_calender_time()
{
    time_t rawtime;
    time (&rawtime);
    char* time_buf_ptr= ctime (&rawtime);
    std::string s(time_buf_ptr);

    std::vector<std::string> tokens;
    std::string token;
    for(int i=0; i<s.size(); i++)
    {
        if(s[i] == ' ' || s[i] == '\n')
        {
            if (!token.empty())
            {
                tokens.push_back(token);
                token.clear();
            }
        }
        else
        {
            token += s[i];
        }
    }

    std::string & day_time = tokens[3];
    int j =0;
    for (int i=0; i<day_time.size(); i++)
    {
            if(day_time[i]== ':' && j == 0)
            {
                    day_time[i]='h';
                    j ++;
            }
            else if (day_time[i]== ':' && j == 1)
            {
                    day_time[i]='m';
                    j ++;
            }
            else if (day_time[i]== ':' && j == 2)
            {
                    day_time[i]='s';
                    break;
            }
    }
    day_time.push_back('s');
    std::string out_str = "date_" + tokens[4]+ '_' +  tokens[1] + '_' +  tokens[2] + '_' +tokens[3] ;
    std::cout << "current date " << out_str << std::endl << std::endl;

    return out_str;
}

/*
void set_input_pair_path(s &input_kmer_name , s &input_kmer_name_2,
                         s &input_read_name,  s &input_read_name_2)
{
    input_kmer_path = base_path_str + "/test_data" + input_kmer_name;
    input_kmer_path_2 = base_path_str + "/test_data" + input_kmer_name_2;
    input_read_path = base_path_str + "/test_data" + input_read_name;
    input_read_path_2 = base_path_str + "/test_data" + input_read_name_2;
    has_pair = true;
}
void set_if_use_pair(bool is_use_pair)
{
    has_pair = is_use_pair;
}
void set_output_dir(s &output_dir_name)
{
    output_path_str = base_path_str + output_dir_name;
    output_contig_path = output_path_str + contig_name;
    output_comp_path = output_path_str + comp_name;
    output_kmer_path = output_path_str + kmer_name;
    deleted_contig_path = output_path_str + "/deleted_kmer";

    output_components_read_dir = output_path_str + "/components_reads";
    output_components_kmer_dir = output_path_str + "/components_kmer";
}
*/
