/* 
 * File:   main.cpp
 * Author: bx
 *
 * Created on August 15, 2017, 2:47 PM
 */
#include "Kmer_handler.h"
#include <unistd.h>
#include "Contig_graph_handler.h"
#include "encoding.h"
#include "local_file_structure.h"
#include "log.h"
#include <stdlib.h>

void test_encoding_decoding();
void test_Kmer_without_filter(Local_files * lf);
void test_Kmer(Local_files * lf);
void test_Contig_graph(Local_files * lf);
void test_load_contig_graph(Local_files * lf);
void test_load_dump(Local_files * lf);
void test_pair_read_contig_graph(Local_files * lf);
void test_load_pair_contig_graph(Local_files * lf);


int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int get_proc_mem_value(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

struct Run_setting {
    Run_setting () 
    {
        kmer_length = 25;
        rmer_length = 15;
        threshold = 0.5;
        min_count = 3;
        min_len = 75;
        is_use_set = false;
    }    
    uint16_t contig_num;
    uint8_t kmer_length;
    uint8_t rmer_length;
    double threshold;
    kmer_count_t min_count; //5
    Contig_handler::size_type min_len; //20
    bool is_use_set;
};

Block_timer timer;



int main(int argc, char** argv) {                
    //char letter;          
    
    int desiredTest = 0;
    while(desiredTest<1 || desiredTest>8)
    {
        std::cout << "Select desired application: \n";
        std::cout << "   1) test_encoding_decoding\n";
        std::cout << "   2) test_Kmer_without_filter\n";
        std::cout << "   3) test_Kmer with list\n";
        std::cout << "   4) test_Contig_graph\n";
        std::cout << "   5) test load and process contig_graph\n";  
        std::cout << "   6) test_load_dump\n"; 
        std::cout << "   7) test pair read contig graph\n";
        std::cout << "   8) test load and process pair read contig graph\n";
        
        std::cin >> desiredTest;
    }
        
    boost::filesystem::path base_path = boost::filesystem::current_path();
    std::string base_path_str(base_path.c_str());            
    
    // file system setting
    std::string output_dir_name("/output");
    std::string input_kmer_name("/medium.dict");     
    std::string input_read_name("/medium.fasta");
    std::string out_log_name("/default_shannonC_log");
    std::string out_contig_name("/contig");
    std::string out_comp_name("/component_array");
    std::string out_kmer_name("/kmer");
        
    Local_files local_files(base_path_str, output_dir_name, input_kmer_name, input_read_name, 
            out_log_name, out_contig_name, out_comp_name, out_kmer_name);                                        
    memcpy(shc_logname, local_files.log_filename_path.c_str(), 
                                        local_files.log_filename_path.size());    
                     
    shc_log_info(shc_logname, "Shannon C is starting\n");
    switch(desiredTest){
        case 1:
            test_encoding_decoding();
            break;
	case 2:
            test_Kmer_without_filter(&local_files);
            break;
        case 3:
            test_Kmer(&local_files);
            break;            
        case 4:
            test_Contig_graph(&local_files);
            break;       
        case 5:
            test_load_contig_graph(&local_files);
            break;
        case 6:            
            test_load_dump(&local_files);
        case 7:            
            test_pair_read_contig_graph(&local_files);
        case 8:
            test_load_pair_contig_graph(&local_files);
            
        default:
            break;
    }                      
    
    //printf("press any key to end\n");
    //std::cin >> letter;
    return 0;
}

void test_load_pair_contig_graph(Local_files * lf)
{
    bool use_multiple_partition = false;
    idx_t partition_size = 500;
    idx_t penalty = 5;
    real_t overload = 1.5;       
    Metis_setup metis_setup(use_multiple_partition, partition_size, penalty, overload);
    
    std::string kmer1 = "/medium_p1.dict", kmer2 = "/medium_p2.dict";
    std::string read1 = "/medium_p1.fasta", read2 = "/medium_p2.fasta";    
    std::string out_dir("/output_pair");
    lf->set_input_pair_path(kmer1, kmer2, read1, read2);
    lf->set_output_dir(out_dir);
    print_local_file_system(lf);
    
    shc_log_info(shc_logname, "Shannon C loads previous files\n");
    Kmer_handler kmer_handler(lf);
    Contig_handler contig_handler; 
    kmer_handler.get_contig_handler(&contig_handler);
    
    kmer_handler.load_kmers_with_info_from_file(lf->output_kmer_path);        
    
    contig_handler.load_contig_file(lf->output_contig_path);
    
    Contig_graph_handler graph_handler(kmer_handler.get_kmer_length()-1, 
                                       &kmer_handler, 
                                       &contig_handler, lf, metis_setup);
    
    graph_handler.group_components();   
    graph_handler.dump_component_array(lf->output_comp_path);
    graph_handler.assign_paired_read_to_components(3);
    graph_handler.assign_kmer_to_components();
}

void test_pair_read_contig_graph(Local_files * lf)
{  
    char letter;
    bool is_list = true;
    std::cout << "If use list for redundancy removal (Y/N)" << std::endl;
    if(std::cin >> letter)
    {
        if(letter=='N' || letter=='n')
        {
            is_list = false;
            std:: cout << "using a set to remove redundancy" << std::endl;
        }
    }    
    // parameter setting
    start_timer(&timer);                    
    std::string kmer1 = "/medium_p1.dict", kmer2 = "/medium_p2.dict";
    std::string read1 = "/medium_p1.fasta", read2 = "/medium_p2.fasta";    
    std::string out_dir("/output_pair");
    lf->set_input_pair_path(kmer1, kmer2, read1, read2);
    lf->set_output_dir(out_dir);
    print_local_file_system(lf);
         
    bool is_compress = false;
    uint8_t kmer_length = 25;
    uint8_t rmer_length = 15;
    double threshold = 0.5;
    kmer_count_t min_count = 3; //5
    Contig_handler::size_type min_len = 75; //20
    bool is_use_set = !is_list;
    
    bool use_multiple_partition = false;
    idx_t partition_size = 500;
    idx_t penalty = 5;
    real_t overload = 1.5;       
    Metis_setup metis_setup(use_multiple_partition, partition_size, penalty, overload);
         
    Kmer_handler kmer_handler(kmer_length,   
            min_count, min_len, threshold, rmer_length, is_use_set, lf);        
    printf("here-1\n");
    Contig_handler contig_handler(is_compress);   
    Contig_graph_handler graph_handler(kmer_length-1, 
                        &kmer_handler, &contig_handler, lf, metis_setup); 
          
    //start processing
    kmer_handler.get_contig_handler(&contig_handler);     
    std::cout << "Start reading kmer file"  << std::endl;
    if(kmer_handler.build_dict_from_kmer_file() != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
    std::cout << "Finish build kmer, with mem " << std::endl;    
    kmer_handler.sort_kmer_descending_count(); 
    std::cout << "after sorting, mem " << std::endl; //get_proc_mem_value()-before_mem << 
    kmer_handler.find_contig();
    std::cout << "after find contig, mem " <<  std::endl;   
    kmer_handler.dump_kmers_to_file(lf->output_kmer_path); 
    contig_handler.dump_all_contig(lf->output_contig_path);    
     
    graph_handler.group_components();   
    
    std::cout << "after group contig " <<  std::endl;
        
    graph_handler.assign_reads_to_components(3);
    graph_handler.assign_kmer_to_components(); 
    
    graph_handler.dump_component_array(lf->output_comp_path);
    
}

void test_load_dump(Local_files * lf)
{
    print_local_file_system(lf);
    std::string kmer_back_file((lf->output_kmer_path) + "_re");
    std::string contig_back_file((lf->output_contig_path) + "_re");
    Kmer_handler kmer_handler(lf);    
    kmer_handler.load_kmers_with_info_from_file(lf->output_kmer_path);
    kmer_handler.dump_kmers_to_file(kmer_back_file);
    
    Contig_handler contig_handler; 
    contig_handler.load_contig_file(lf->output_contig_path);
    contig_handler.dump_all_contig(contig_back_file);    
}

void test_load_contig_graph(Local_files * lf)
{
    bool use_multiple_partition = false;
    idx_t partition_size = 500;
    idx_t penalty = 5;
    real_t overload = 1.5;       
    Metis_setup metis_setup(use_multiple_partition, partition_size, penalty, overload);
    
    print_local_file_system(lf);
    shc_log_info(shc_logname, "Shannon C loads previous files\n");
    Kmer_handler kmer_handler(lf);
    Contig_handler contig_handler; 
    kmer_handler.get_contig_handler(&contig_handler);
    
    kmer_handler.load_kmers_with_info_from_file(lf->output_kmer_path);        
    
    contig_handler.load_contig_file(lf->output_contig_path);
    
    Contig_graph_handler graph_handler(kmer_handler.get_kmer_length()-1, 
                                       &kmer_handler, 
                                       &contig_handler, lf, metis_setup);
    
    graph_handler.group_components();   
    graph_handler.dump_component_array(lf->output_comp_path);
    graph_handler.assign_reads_to_components(3);
    graph_handler.assign_kmer_to_components();
    
}

void test_Contig_graph(Local_files * lf)
{  
    char letter;
    bool is_list = true;
    bool is_compress = false;
    print_local_file_system(lf);
    std::cout << "If use list for redundancy removal (Y/N)" << std::endl;
    if(std::cin >> letter)
    {
        if(letter=='N' || letter=='n')
        {
            is_list = false;
            std:: cout << "using a set to remove redundancy" << std::endl;
        }
    }
    
    // parameter setting
    start_timer(&timer);
    
    uint16_t contig_num;
    uint8_t kmer_length = 25;
    uint8_t rmer_length = 15; 
    double threshold = 0.5;
    kmer_count_t min_count = 3; //5
    Contig_handler::size_type min_len = 75; //20
    bool is_use_set = !is_list;
    //metis_setup
    bool use_multiple_partition = false;
    idx_t partition_size = 500;
    idx_t penalty = 5;
    real_t overload = 1.5;       
    Metis_setup metis_setup(use_multiple_partition, partition_size, penalty, overload);
    
    
    std::cout << "using set? " << (is_use_set? "Yes":"no") << std::endl;
    //int before_mem = get_proc_mem_value();
    
    Kmer_handler kmer_handler(kmer_length,   
            min_count, min_len, threshold, rmer_length, is_use_set, lf);            
    
    Contig_handler contig_handler(is_compress);       
    kmer_handler.get_contig_handler(&contig_handler); 
    
    std::cout << "Start reading kmer file"  << std::endl;
    if(kmer_handler.build_dict_from_kmer_file() != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
    std::cout << "Finish build kmer, with mem " << std::endl;    
    kmer_handler.sort_kmer_descending_count(); 
    std::cout << "after sorting, mem " << std::endl; //get_proc_mem_value()-before_mem << 
    contig_num = kmer_handler.find_contig();
    std::cout << "after find contig, mem " <<  std::endl;
    
    kmer_handler.dump_kmers_to_file(lf->output_kmer_path);
    contig_handler.dump_all_contig(lf->output_contig_path);
    
    Contig_graph_handler graph_handler(kmer_length-1, 
                        &kmer_handler, &contig_handler, lf, metis_setup);
    graph_handler.group_components();   
    
    std::cout << "after group contig " <<  std::endl;
    
    graph_handler.dump_component_array(lf->output_comp_path);
    graph_handler.assign_reads_to_components(3);
    graph_handler.assign_kmer_to_components();    
    
    std::cout << "The whole process finish ";
    stop_timer(&timer);
}

void test_Kmer(Local_files * lf)
{  
    char letter;
    bool is_list = true;
    std::cout << "If use list for redundancy removal (Y/N)" << std::endl;
    if(std::cin >> letter)
    {
        if(letter=='N' || letter=='n')
        {
            is_list = false;
            std:: cout << "using a set to remove redundancy" << std::endl;
        }
    }
    print_local_file_system(lf);
    uint64_t contig_num;   
    uint8_t kmer_length = 25;
    uint8_t rmer_length = 15;
    double threshold = 0.5; //0.5
    kmer_count_t min_count = 1; //5
    Contig_handler::size_type min_len = 75; //20
    bool is_use_set = !is_list;        
    
    Kmer_handler kmer_handler(kmer_length,  
            min_count, min_len, threshold, rmer_length, is_use_set, lf);     
    
    Contig_handler contig_handler;       
    kmer_handler.get_contig_handler(&contig_handler);    
    
    if(kmer_handler.build_dict_from_kmer_file() != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
    kmer_handler.sort_kmer_descending_count(); 
    contig_num = kmer_handler.find_contig();
    printf("number of contig is %" PRIu64 "\n", contig_num);  
    //kmer_handler.dump_kmers_to_file(lf->kmer_output_path);    
    //contig_handler.dump_all_contig(lf->contig_output_path);
}



void test_Kmer_without_filter(Local_files * lf)
{    
    char letter;
    bool is_list = true;
    std::cout << "If use list for redundancy removal (Y/N)" << std::endl;
    if(std::cin >> letter)
    {
        if(letter=='N' || letter=='n')
        {
            is_list = false;
            std:: cout << "using a set to remove redundancy" << std::endl;
        }
    }

    print_local_file_system(lf);
    uint64_t contig_num;
    
    uint8_t kmer_length = 25;
    uint8_t rmer_length = 15;
    double threshold = 1;
    kmer_count_t min_count = 0; //5
    Contig_handler::size_type min_len = 0; //20
    bool is_use_set = !is_list;
    

    Kmer_handler kmer_handler(kmer_length, 
            min_count, min_len, threshold, rmer_length, is_use_set, lf);     
    
    Contig_handler contig_handler;       
    kmer_handler.get_contig_handler(&contig_handler);
     
    if(kmer_handler.build_dict_from_kmer_file() != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
    kmer_handler.sort_kmer_descending_count(); 
    contig_num = kmer_handler.find_contig();
    printf("number of contig is %" PRIu64 "\n", contig_num);    
    contig_handler.print_delimitor();
    for(size_t i=0; i<contig_handler.delimitor[1];i++)
    {
        printf("%u ", contig_handler.contig_list[i]);
    }        
    //contig_handler.dump_all_contig(output_contig_filename);
}


/**
 * 
 */
void test_encoding_decoding()
{
    uint8_t byte_list[4];
    const char base_string[20] = "CCCCAAAATTTTGGGG";
    char back[20];
    uint16_t length = 16;
    back[length] = '\0';
    encode_base_string(base_string, byte_list , length);
    decode_byte_list(byte_list , back, length);
    if(strcmp(back, base_string) == 0 )    
        printf("CORRECT, encode_base_string decode_base_string works \n");
    else
        printf("WRONG, encode_base_string decode_base_string\n");
        
    uint64_t combined = combine_byte(byte_list, 4);
                
    uint64_t byte = 0;
    encode_kmer(base_string, &byte, length);
        
    if(combined == byte)
        printf("CORRECT, combining bytes equals to kmer direct encode\n");
    else
        printf("WRONG, combining bytes not equals to kmer direct encode\n");
    
    encode_kmer( base_string, &byte, length);
    decode_kmer( back, &byte, length);
    if(strcmp(back, base_string) == 0 )    
        printf("CORRECT, encode/ decoder kmer directly works\n");
    else
        printf("WRONG, encode_base_string decode_base_string\n");
    
}
