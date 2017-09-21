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

using namespace std;

void test_encoding_decoding();
void test_Kmer_without_filter();
void test_Kmer();
void test_Contig_graph();
void test_all();
void test_xmer_at_byte();

char shc_logname[100];

int main(int argc, char** argv) {                
    char letter;         
    
    int desiredTest = 4;
    
    switch(desiredTest){
        case 1:
            test_encoding_decoding();
            break;
	case 2:
            test_Kmer_without_filter();
            break;
        case 3:
            test_Kmer();
            break;            
        case 4:
            test_Contig_graph();
            break;
        case 5:
            test_all();
            break;            
        case 6:
            test_xmer_at_byte();
            break;   
        default:
            break;
    }                      
    
    printf("press any key to end\n");
    std::cin >> letter;
    return 0;
}

void test_xmer_at_byte()
{
    uint8_t xmer_length = 5;
    uint8_t byte_list[3];
    char base[11] = "ACTTCAGGAT";
    char back[11];
    back[10]='\0';    
    encode_base_string( base, byte_list, 10);
    decode_byte_list( byte_list, back , 10);
    printf("%s\n",back);
    
    back[xmer_length]='\0';    
    for(int i=0; i<6;i++)
    {
        get_xmer_at_index(byte_list,xmer_length, i, xmer_length, back);    
        printf("%d: %s\n", i, back);
    }        
}
    
void test_all()
{
    /*
    char cwd[1024];
    uint16_t contig_num;
    uint8_t kmer_length = 25;
    uint8_t rmer_length = 12;
    double threshold = 0.5;
    kmer_count_t min_count = 5;
    Contig_handler::size_type min_len = 20;
    
    
    if (getcwd(cwd, sizeof(cwd)) != NULL)
        printf("Current work dir: %s\n", cwd);
    else
        printf("getcwd error");
        
    std::string input_filename("./test_data/k1mer.dict");
    std::string output_contig_filename("./output/filter_contig.txt");
    
    shc_log_info(shc_logname,"Shannon C is starting\n");
    Kmer_handler kmer_handler(kmer_length, shc_logname);     
    
    Contig_handler contig_handler(shc_logname, rmer_length);       
    kmer_handler.get_contig_handler(&contig_handler);
     
    if(kmer_handler.read_sequence_from_file(input_filename) != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
    kmer_handler.sort_kmer_descending_count(); 
    contig_num = kmer_handler.find_contig();
    
    //contig_handler.filter_contig(min_count, min_len, threshold);
        
    Contig_graph_handler graph_handler(contig_num, kmer_length-1, &kmer_handler,
                                             &contig_handler, shc_logname);
    graph_handler.construct_graph();    
    //graph_handler.record_and_delete_isolated_contig();
        
    //contig_handler.dump_all_contig(output_contig_filename);
    */
}

void test_Contig_graph()
{    
    char cwd[1024];
    uint16_t contig_num;
    uint8_t kmer_length = 25;
        
    typedef boost::filesystem::path path_t;
    path_t base_path = boost::filesystem::current_path();
    std::string base_path_str(base_path.c_str());
    std::string input_kmer_filename = base_path_str + "/test_data/medium_kmer.dict";
    std::string input_read_filename = base_path_str + "/test_data/reads.fasta";
    std::string log_filename_path = base_path_str + "/log/default_shannonC_log";
    std::string contig_output_path = base_path_str + "/output/contig";
    std::string comp_output_path = base_path_str + "/output/component_array";
    
    memcpy(shc_logname, log_filename_path.c_str(), log_filename_path.size());    
    
    std::cout << "current work dir is " << base_path << std::endl<< std::endl;
        
    uint8_t rmer_length = 15;
    double threshold = 0.5;
    kmer_count_t min_count = 1; //5
    Contig_handler::size_type min_len = 75; //20
    bool is_use_set = false;
                
    shc_log_info(shc_logname, "Shannon C is starting\n");
    Kmer_handler kmer_handler(kmer_length, shc_logname,  min_count, 
                min_len, threshold, 
                rmer_length, is_use_set);     
    
    Contig_handler contig_handler(shc_logname);       
    kmer_handler.get_contig_handler(&contig_handler);    
    
    if(kmer_handler.read_sequence_from_file(input_kmer_filename) != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
    kmer_handler.sort_kmer_descending_count(); 
    contig_num = kmer_handler.find_contig();
    //kmer_handler.traverse_kmer_count();
    //std::string contig_output_filename = output_dir + "contig.txt";
    //contig_handler.dump_all_contig(contig_output_path);
    
    Contig_graph_handler graph_handler(kmer_length-1, &kmer_handler,
                             &contig_handler, shc_logname);
    graph_handler.group_components();   
    graph_handler.dump_component_array(comp_output_path);
    graph_handler.assign_reads_to_components(input_read_filename, 3);
    graph_handler.assign_kmer_to_components();
    //graph_handler.record_and_delete_isolated_contig();
    //graph_handler.create_metis_format_from_graph();
    //graph_handler.print_adjacency_structure();
    //graph_handler.add_weighted_edge();
}

void test_Kmer()
{    
    uint64_t contig_num;   
    uint8_t kmer_length = 25;
    uint8_t rmer_length = 15;
    double threshold = 0.5; //0.5
    kmer_count_t min_count = 1; //5
    Contig_handler::size_type min_len = 75; //20
    bool is_use_set = false;
        
    typedef boost::filesystem::path path_t;
    path_t base_path = boost::filesystem::current_path();
    std::string base_path_str(base_path.c_str());
    std::string input_kmer_filename = base_path_str + "/test_data/kmer_read.dict";
    std::string input_read_filename = base_path_str + "/test_data/SE_read.fasta";
    std::string log_filename_path = base_path_str + "/log/default_shannonC_log";
    std::string contig_output_path = base_path_str + "/output/contig";
    
    memcpy(shc_logname, log_filename_path.c_str(), log_filename_path.size());    
    
    std::cout << "current work dir is " << base_path << std::endl<< std::endl;        
    
    shc_log_info(shc_logname,"Shannon C is starting\n");
    Kmer_handler kmer_handler(kmer_length, shc_logname,  min_count, 
                min_len, threshold, 
                rmer_length, is_use_set);     
    
    Contig_handler contig_handler(shc_logname);       
    kmer_handler.get_contig_handler(&contig_handler);
    //kmer_handler.estimate_num_kmer(input_kmer_filename);
    
    if(kmer_handler.read_sequence_from_file(input_kmer_filename) != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
    kmer_handler.sort_kmer_descending_count(); 
    contig_num = kmer_handler.find_contig();
    printf("number of contig is %" PRIu64 "\n", contig_num);  
    
    //contig_handler.print_delimitor();
    //for(size_t i=0; i<contig_handler.delimitor[1];i++)
    //{
    //   printf("%u ", contig_handler.contig_list[i]);
    //}
        
    contig_handler.dump_all_contig(contig_output_path);
}



void test_Kmer_without_filter()
{    
    uint64_t contig_num;
    
    uint8_t kmer_length = 25;
    uint8_t rmer_length = 15;
    double threshold = 1;
    kmer_count_t min_count = 0; //5
    Contig_handler::size_type min_len = 0; //20
    bool is_use_set = false;
    
    typedef boost::filesystem::path path_t;
    path_t base_path = boost::filesystem::current_path();
    std::string base_path_str(base_path.c_str());
    std::string input_kmer_filename = base_path_str + "/test_data/kmer_read.dict";
    std::string input_read_filename = base_path_str + "/test_data/SE_read.fasta";
    std::string log_filename_path = base_path_str + "/log/default_shannonC_log";
    
    memcpy(shc_logname, log_filename_path.c_str(), log_filename_path.size());           
    
    std::cout << "current work dir is " << base_path << std::endl<< std::endl;
            
    shc_log_info(shc_logname,"Shannon C is starting\n");
    Kmer_handler kmer_handler(kmer_length, shc_logname,  min_count, 
                min_len, threshold, 
                rmer_length, is_use_set);     
    
    Contig_handler contig_handler(shc_logname);       
    kmer_handler.get_contig_handler(&contig_handler);
     
    if(kmer_handler.read_sequence_from_file(input_kmer_filename) != 0)
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

/* old tests */
/*
 char base[5];

    base[0] = 'T';
    base[1] = 'T';
    base[2] = 'A';
    base[3] = 'A';
    base[4] = '\0';
  
    
    char file_name[50] = "base_string.txt"; 
    char base_string[31];
    char back_string[31];
    back_string[30] = '\0';
    uint8_t byte_list[10] = {0};
    read_input(base_string, file_name);        
    
    uint16_t out_length;
    double in_length = 30;
    out_length = encode_base_string(base_string, byte_list, in_length);
    
    printf("out length %u\n", out_length);
    print_byte_list(byte_list, out_length);
    decode_byte_list(byte_list, back_string, in_length);
    printf("%s\n", back_string);
 
 */
