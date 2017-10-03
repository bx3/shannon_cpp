#include "local_file_structure.h"

void replace_directory(boost::filesystem::path & dir_path)
{
    boost::filesystem::remove_all(dir_path);
    add_directory(dir_path);
}

void add_directory(boost::filesystem::path & dir_path)
{     
    if(boost::filesystem::create_directory(dir_path))
    {
        std::cout << "success to add directory: " << dir_path << std::endl;
    }
    else
    {
        std::cout << "fail to add directory: " << dir_path << std::endl;
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

void print_local_file_system(Local_files *local_files)
{
    std::cout << "\033[1;32m";  //34 blue, 31 red, 35 purple 
    std::cout << "current work dir is      : " << local_files->base_path_str << std::endl;
    
    std::cout << "Input                     *********** " << std::endl;
    std::cout << "If use pair end read     : " << ((local_files->has_pair)?"Yes":"No") << std::endl;
    std::cout << "kmer input from          : " << local_files->input_kmer_path << std::endl;
    if(local_files->has_pair)
        std::cout << "kmer pair input from     : " << local_files->input_kmer_path_2 << std::endl;    
    std::cout << "read input from          : " << local_files->input_read_path << std::endl;
    if(local_files->has_pair)
        std::cout << "read input from          : " << local_files->input_read_path_2 << std::endl;
    
    std::cout << "Output                    *********** " << std::endl;
    std::cout << "kmer output save to      : " << local_files->output_kmer_path << std::endl;
    std::cout << "contig output saved to   : " << local_files->output_contig_path << std::endl;
    std::cout << "component output saved to: " << local_files->output_comp_path << std::endl;
    
    std::cout << "Others                    *********** " << std::endl;
    std::cout << "log file save to         : " << local_files->log_filename_path << std::endl;
    std::cout << "deleted kmer saved to    : " << local_files->deleted_contig_path << std::endl;
    std::cout << "\033[0m" << std::endl;   
}