#include "local_file_structure.h"

void add_directory(boost::filesystem::path & dir_path)
{ 
    std::cout << "Trying to add dir: " << dir_path <<std::endl;
    if(!boost::filesystem::exists(dir_path))
    {
        if(boost::filesystem::create_directory(dir_path))
        {
            std::cout << "success to add directory" <<std::endl;
        }
        else
        {
            std::cout << "fail to add directory" << std::endl;
        }
    }
    else
    {
        std::cout << "directory already exists" <<std::endl;
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