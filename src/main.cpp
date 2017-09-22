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

using namespace std;

void test_encoding_decoding();
void test_Kmer_without_filter(Local_files * lf);
void test_Kmer(Local_files * lf);
void test_Contig_graph(Local_files * lf);
void test_load_contig_graph(Local_files * lf);
void test_load_dump(Local_files * lf);


int main(int argc, char** argv) {                
    char letter;             
    
    int desiredTest = 0;
    while(desiredTest<1 || desiredTest>6)
    {
        std::cout << "Select desired application: \n";
        std::cout << "   1) test_encoding_decoding\n";
        std::cout << "   2) test_Kmer_without_filter\n";
        std::cout << "   3) test_Kmer\n";
        std::cout << "   4) test_Contig_graph\n";
        std::cout << "   5) test_load_contig_graph\n";  
        std::cout << "   6) test_load_dump\n"; 
        
        std::cin >> desiredTest;
    }
    
    boost::filesystem::path base_path = boost::filesystem::current_path();
    std::string base_path_str(base_path.c_str());        
    
    // file system setting
    std::string input_kmer_name("/medium.dict");     
    std::string input_read_name("/medium.fasta");
    std::string out_log_name("/default_shannonC_log");
    std::string out_contig_name("/contig");
    std::string out_comp_name("/component_array");
    std::string out_kmer_name("/kmer");
        
    Local_files local_files(base_path_str, input_kmer_name, input_read_name, 
            out_log_name, out_contig_name, out_comp_name, out_kmer_name);                                        
    memcpy(shc_logname, local_files.log_filename_path.c_str(), 
                                        local_files.log_filename_path.size());
    std::cout << std::endl;
    std::cout << "\033[1;32m";  //34 blue, 31 red, 35 purple 
    std::cout << "current work dir is      : " << local_files.base_path_str << std::endl;
    std::cout << "kmer input from          : " << local_files.input_kmer_path << std::endl;
    std::cout << "read input from          : " << local_files.input_read_path << std::endl;
    std::cout << "log file save to         : " << local_files.log_filename_path << std::endl;
    std::cout << "kmer output save to      : " << local_files.kmer_output_path << std::endl;
    std::cout << "contig output saved to   : " << local_files.contig_output_path << std::endl;
    std::cout << "component output saved to: " << local_files.comp_output_path << std::endl << std::endl;
    std::cout << "\033[0m";
         
    
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
        default:
            break;
    }                      
    
    printf("press any key to end\n");
    std::cin >> letter;
    return 0;
}

void test_load_dump(Local_files * lf)
{
    std::string kmer_back_file((lf->kmer_output_path) + "_re");
    std::string contig_back_file((lf->contig_output_path) + "_re");
    Kmer_handler kmer_handler;    
    kmer_handler.load_kmers_with_info_from_file(lf->kmer_output_path);
    kmer_handler.dump_kmers_to_file(kmer_back_file);
    
    Contig_handler contig_handler; 
    contig_handler.load_contig_file(lf->contig_output_path);
    contig_handler.dump_all_contig(contig_back_file);
    
}

void test_load_contig_graph(Local_files * lf)
{
    shc_log_info(shc_logname, "Shannon C loads previous files\n");
    Kmer_handler kmer_handler;
    Contig_handler contig_handler; 
    kmer_handler.get_contig_handler(&contig_handler);
    
    kmer_handler.load_kmers_with_info_from_file(lf->kmer_output_path);
    
    std::cout << "inside main " << (uint16_t)kmer_handler.get_kmer_length() << std::endl;
    
    contig_handler.load_contig_file(lf->contig_output_path);
    
    Contig_graph_handler graph_handler(kmer_handler.get_kmer_length()-1, 
                                       &kmer_handler, 
                                       &contig_handler);
    
    graph_handler.group_components();   
    graph_handler.dump_component_array(lf->comp_output_path);
    graph_handler.assign_reads_to_components(lf->input_read_path, 3);
    graph_handler.assign_kmer_to_components();
    
}

void test_Contig_graph(Local_files * lf)
{  
    // parameter setting
    uint16_t contig_num;
    uint8_t kmer_length = 25;
    uint8_t rmer_length = 15;
    double threshold = 0.5;
    kmer_count_t min_count = 1; //5
    Contig_handler::size_type min_len = 75; //20
    bool is_use_set = false;
    
    Kmer_handler kmer_handler(kmer_length,   
            min_count, min_len, threshold, rmer_length, is_use_set);     
    
    Contig_handler contig_handler;       
    kmer_handler.get_contig_handler(&contig_handler);    
    
    if(kmer_handler.read_sequence_from_file(lf->input_kmer_path) != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
    kmer_handler.sort_kmer_descending_count(); 
    contig_num = kmer_handler.find_contig();
    
    kmer_handler.dump_kmers_to_file(lf->kmer_output_path);
    contig_handler.dump_all_contig(lf->contig_output_path);
    
    Contig_graph_handler graph_handler(kmer_length-1, &kmer_handler, &contig_handler);
    graph_handler.group_components();   
    graph_handler.dump_component_array(lf->comp_output_path);
    graph_handler.assign_reads_to_components(lf->input_read_path, 3);
    graph_handler.assign_kmer_to_components();    
}

void test_Kmer(Local_files * lf)
{    
    uint64_t contig_num;   
    uint8_t kmer_length = 25;
    uint8_t rmer_length = 15;
    double threshold = 0.5; //0.5
    kmer_count_t min_count = 1; //5
    Contig_handler::size_type min_len = 75; //20
    bool is_use_set = false;        
    
    Kmer_handler kmer_handler(kmer_length,  
            min_count, min_len, threshold, rmer_length, is_use_set);     
    
    Contig_handler contig_handler;       
    kmer_handler.get_contig_handler(&contig_handler);    
    
    if(kmer_handler.read_sequence_from_file(lf->input_kmer_path) != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
    kmer_handler.sort_kmer_descending_count(); 
    contig_num = kmer_handler.find_contig();
    printf("number of contig is %" PRIu64 "\n", contig_num);  
        
    contig_handler.dump_all_contig(lf->contig_output_path);
}



void test_Kmer_without_filter(Local_files * lf)
{    
    uint64_t contig_num;
    
    uint8_t kmer_length = 25;
    uint8_t rmer_length = 15;
    double threshold = 1;
    kmer_count_t min_count = 0; //5
    Contig_handler::size_type min_len = 0; //20
    bool is_use_set = false;
    

    Kmer_handler kmer_handler(kmer_length, 
            min_count, min_len, threshold, rmer_length, is_use_set);     
    
    Contig_handler contig_handler;       
    kmer_handler.get_contig_handler(&contig_handler);
     
    if(kmer_handler.read_sequence_from_file(lf->input_kmer_path) != 0)
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
