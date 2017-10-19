#include "json_parser.h"

void parser_setting_file(std::string & file_path, Shannon_C_setting & setting)
{
    using boost::property_tree::ptree;    

    std::ifstream json_file(file_path.c_str());

    ptree pt;

    boost::property_tree::json_parser::read_json(json_file, pt);
    // top level
    uint8_t kmer_length = pt.get<uint8_t>("kmer_length");
    bool has_pair = pt.get<bool>("has_pair");
    bool is_compress = pt.get<bool>("is_compress");
    
    setting.kmer_length = kmer_length;
    setting.has_pair = has_pair;
    setting.is_compress = is_compress;
                
    // file system
    ptree file_structure = pt.get_child("file_structure");
    std::string input_kmer_name = file_structure.get<std::string>("input_kmer_name");
    std::string input_kmer_name_2 = file_structure.get<std::string>("input_kmer_name_2");
    std::string input_read_name = file_structure.get<std::string>("input_read_name");
    std::string input_read_name_2 = file_structure.get<std::string>("input_read_name_2");
    std::string output_dir_name = file_structure.get<std::string>("output_dir_name");
    std::string out_log_name = file_structure.get<std::string>("out_log_name");
    std::string out_contig_name = file_structure.get<std::string>("out_contig_name");
    std::string out_comp_name = file_structure.get<std::string>("out_comp_name");
    std::string out_kmer_name = file_structure.get<std::string>("out_kmer_name");    
    
    boost::filesystem::path base_path = boost::filesystem::current_path();
    std::string base_path_str(base_path.c_str());  
    
    if(input_kmer_name_2.empty() && input_read_name_2.empty())
    {
        if(has_pair)
        {
            shc_log_error("has_pair setting is inconsistent to file input");
        }
        
        Local_files local_files(
             base_path_str, output_dir_name, input_kmer_name, input_read_name, 
            out_log_name, out_contig_name, out_comp_name, out_kmer_name); 
        setting.local_files = local_files;
    }
    else
    {
        if(!has_pair)
        {
            shc_log_error("has_pair setting is inconsistent to file input");
        }
                
        Local_files local_files(
                base_path_str, output_dir_name, input_kmer_name, 
                input_kmer_name_2, input_read_name, input_read_name_2,
                out_log_name, out_contig_name, out_comp_name, out_kmer_name);         
        setting.local_files = local_files;
    }
    
    
        
    // dup_correction
    ptree dup_correction = pt.get_child("dup_correction");
    uint8_t rmer_length = dup_correction.get<uint8_t>("rmer_length");
    double threshold = dup_correction.get<double>("threshold");
    kmer_count_t min_count = dup_correction.get<kmer_count_t>("min_count");
    size_t min_length = dup_correction.get<size_t>("min_length");
    bool is_use_set = dup_correction.get<bool>("is_use_set");
    
    setting.dup_setting.set_para(rmer_length, threshold, min_count, min_length, is_use_set);    
            
    //metis_setup
    ptree metis_setup = pt.get_child("metis_setup");
    bool use_multiple_partition = metis_setup.get<bool>
                                                    ("use_multiple_partition");
    idx_t partition_size = metis_setup.get<idx_t>("partition_size");
    idx_t penalty = metis_setup.get<idx_t>("penalty");
    real_t overload = metis_setup.get<real_t>("overload");
        
    setting.metis_setup.set_para(use_multiple_partition, partition_size, penalty, 
                                                                overload);  
    
    ptree output_setup = pt.get_child("output_setup");
    bool kmer_with_info = output_setup.get<bool>("kmer_with_info");
    bool contig = output_setup.get<bool>("contig");
    bool component_array = output_setup.get<bool>("component_array");
    bool read_to_component = output_setup.get<bool>("read_to_component");
    bool kmer_to_component = output_setup.get<bool>("kmer_to_component");
    
    setting.output_setup.set_para(kmer_with_info, contig, component_array, 
                         read_to_component, kmer_to_component);
    
}

void print_setting(Shannon_C_setting & setting)
{        
    std::cout << "\033[1;31m";  //34 blue, 31 red, 35 purple       
    std::cout << "General setting" << std::endl;
    std::cout << "kmer_length is           : " << static_cast<uint16_t>(setting.kmer_length) << std::endl;
    std::cout << "has_pair                 : " << ((setting.has_pair)?("Yes"):("No")) << std::endl;
    std::cout << "use compressed contig    : " << ((setting.is_compress)?("Yes"):("No")) << std::endl;
    std::cout << std::endl;
    
    print_local_file_system(&setting.local_files);
    
    std::cout << "\033[1;33m";  //34 blue, 31 red, 35 purple       
    std::cout << "Duplicate correction      *********** " << std::endl;
    std::cout << "rmer_length              : " << (uint16_t)setting.dup_setting.rmer_length << std::endl;
    std::cout << "min_count                : " << setting.dup_setting.min_count << std::endl;
    std::cout << "min_len                  : " << setting.dup_setting.min_len << std::endl;
    std::cout << "is_use_set               : " << ((setting.dup_setting.is_use_set)?("Yes"):("No")) << std::endl;
    std::cout << std::endl;
    
    std::cout << "\033[1;35m";  //34 blue, 31 red, 35 purple       
    std::cout << "Metis setup              *********** " << std::endl;
    std::cout << "use_multiple_partition   : " << ((setting.metis_setup.is_multiple_partition)?("Yes"):("No")) << std::endl;    
    std::cout << "partition_size           : " << setting.metis_setup.partition_size << std::endl;
    std::cout << "penalty `                : " << setting.metis_setup.penalty << std::endl;
    std::cout << "overload                 : " << setting.metis_setup.overload << std::endl;       
    std::cout << std::endl;
    
    std::cout << "\033[1;34m";  //34 blue, 31 red, 35 purple       
    std::cout << "Output setting           *********** " << std::endl;
    std::cout << "dump kmer with info      : " << ((setting.output_setup.kmer_with_info)?("Yes"):("No")) << std::endl;    
    std::cout << "dump contig              : " << ((setting.output_setup.contig)?("Yes"):("No")) << std::endl;    
    std::cout << "dump component_array     : " << ((setting.output_setup.component_array)?("Yes"):("No")) << std::endl;    
    std::cout << "dump read to component   : " << ((setting.output_setup.read_to_component)?("Yes"):("No")) << std::endl;    
    std::cout << "dump kmer to component   : " << ((setting.output_setup.kmer_to_component)?("Yes"):("No")) << std::endl;    
    
    std::cout << "\033[0m" << std::endl;   
    
    
}