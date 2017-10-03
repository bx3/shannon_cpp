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
#include "json_parser.h"

void test_encoding_decoding();
void test_Kmer(Shannon_C_setting * setting);
void test_Contig_graph(Shannon_C_setting * setting);
void test_load_contig_graph(Shannon_C_setting * setting);
void test_load_dump(Shannon_C_setting * setting);
void test_pair_read_contig_graph(Shannon_C_setting * setting);
void test_load_pair_contig_graph(Shannon_C_setting * setting);


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

Block_timer timer;



int main(int argc, char** argv) {                
    
    if(argc != 2)
    {
        shc_log_error("usage: ./Shannon_C_seq setting_file_name\n");
        exit(1);
    }
    
    std::string setting_file_name(argv[1]);
    
    boost::filesystem::path base_path = boost::filesystem::current_path();
    std::string base_path_str(base_path.c_str());            
    //("/shannon_C_setting.txt");    
    std::string setting_path = base_path_str + "/" + setting_file_name;
    
    if(!boost::filesystem::exists(setting_path))
    {
        shc_log_error("No such file: %s\n", setting_path.c_str());
        exit(1);
    }
    
    Shannon_C_setting setting;
        
    parser_setting_file(setting_path, setting);
    print_setting(setting);
    
    
    int desiredTest = 0;
    while(desiredTest<1 || desiredTest>7)
    {
        std::cout << "Select desired application: \n";
        std::cout << "   1) test_encoding_decoding\n";
        std::cout << "   2) test_Kmer with list\n";
        std::cout << "   3) test_Contig_graph\n";
        std::cout << "   4) test load and process contig_graph\n";  
        std::cout << "   5) test_load_dump\n"; 
        std::cout << "   6) test pair read contig graph\n";
        std::cout << "   7) test load and process pair read contig graph\n";
        
        std::cin >> desiredTest;
    }
                                                           
    //memcpy(shc_logname, local_files.log_filename_path.c_str(), 
    //                                    local_files.log_filename_path.size());    
                     
    shc_log_info(shc_logname, "Shannon C is starting\n");
    switch(desiredTest){
        case 1:
            test_encoding_decoding();
            break;
        case 2:
            test_Kmer(&setting);
            break;            
        case 3:
            test_Contig_graph(&setting);
            break;       
        case 4:
            test_load_contig_graph(&setting);
            break;
        case 5:            
            test_load_dump(&setting);
        case 6:            
            test_pair_read_contig_graph(&setting);
        case 7:
            test_load_pair_contig_graph(&setting);
            
        default:
            break;
    }                          
    return 0;
}

void test_load_pair_contig_graph(Shannon_C_setting * setting)
{    
    shc_log_info(shc_logname, "Shannon C loads previous files\n");
    Kmer_handler kmer_handler(&setting->local_files);
    Contig_handler contig_handler; 
    kmer_handler.get_contig_handler(&contig_handler);
    
    kmer_handler.load_kmers_with_info_from_file(setting->local_files.output_kmer_path);        
    
    contig_handler.load_contig_file(setting->local_files.output_contig_path);
    
    Contig_graph_handler graph_handler(kmer_handler.get_kmer_length()-1, 
                             &kmer_handler, &contig_handler,
                              &setting->local_files, setting->metis_setup);
    
    graph_handler.group_components();   
    
    
    if(setting->output_setup.component_array)
        graph_handler.dump_component_array(setting->local_files.output_comp_path);
    
    if(setting->output_setup.read_to_component)
        graph_handler.assign_paired_read_to_components(3);
    
    if(setting->output_setup.kmer_to_component)
        graph_handler.assign_kmer_to_components();
}

void test_pair_read_contig_graph(Shannon_C_setting * setting)
{      
    Duplicate_setting & dup = setting->dup_setting;
    // parameter setting
    start_timer(&timer);                                              
    Kmer_handler kmer_handler(setting->kmer_length,   
            dup.min_count, dup.min_len, dup.threshold, dup.rmer_length, 
            dup.is_use_set, &setting->local_files);        
    
    Contig_handler contig_handler(setting->is_compress);   
    Contig_graph_handler graph_handler(setting->kmer_length-1, &kmer_handler,
                  &contig_handler, &setting->local_files, setting->metis_setup); 
          
    //start processing
    kmer_handler.get_contig_handler(&contig_handler);     
    std::cout << "Start reading kmer file"  << std::endl;
    if(kmer_handler.build_dict_from_kmer_file() != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
      
    kmer_handler.sort_kmer_descending_count(); 
    
    kmer_handler.find_contig();
    
    kmer_handler.dump_kmers_to_file(setting->local_files.output_kmer_path); 
    contig_handler.dump_all_contig(setting->local_files.output_contig_path);    
     
    graph_handler.group_components();   
    
    std::cout << "after group contig " <<  std::endl;
        
    
    if(setting->output_setup.read_to_component)
        graph_handler.assign_reads_to_components(3);
    if(setting->output_setup.kmer_to_component)
        graph_handler.assign_kmer_to_components(); 
    if(setting->output_setup.component_array)
        graph_handler.dump_component_array(setting->local_files.output_comp_path);
    
}

void test_load_dump(Shannon_C_setting * setting)
{
    Local_files * lf = &setting->local_files;
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

void test_load_contig_graph(Shannon_C_setting * setting)
{
    Duplicate_setting & dup = setting->dup_setting;
    Metis_setup & metis_setup = setting->metis_setup;    
    Local_files * lf = &setting->local_files;
    
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
    
    if(setting->output_setup.read_to_component)
        graph_handler.assign_reads_to_components(3);
    
    if(setting->output_setup.read_to_component)
        graph_handler.assign_reads_to_components();    
}

void test_Contig_graph(Shannon_C_setting * setting)
{  
    Duplicate_setting & dup = setting->dup_setting;
    Metis_setup & metis_setup = setting->metis_setup;    
    Local_files * lf = &setting->local_files;
    
    // parameter setting
    start_timer(&timer);            
    
    Kmer_handler kmer_handler(setting->kmer_length, dup.min_count, dup.min_len, 
            dup.threshold, dup.rmer_length, dup.is_use_set, lf);            
    
    Contig_handler contig_handler(setting->is_compress);       
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
    
    Contig_graph_handler graph_handler(setting->kmer_length-1, 
                        &kmer_handler, &contig_handler, lf, metis_setup);
    graph_handler.group_components();   
    
    std::cout << "after group contig " <<  std::endl;
    
    graph_handler.dump_component_array(lf->output_comp_path);
    graph_handler.assign_reads_to_components(3);
    graph_handler.assign_kmer_to_components();    
    
    std::cout << "The whole process finish ";
    stop_timer(&timer);
}

void test_Kmer(Shannon_C_setting * setting)
{  
    Duplicate_setting & dup = setting->dup_setting;
    Metis_setup & metis_setup= setting->metis_setup;    
    Local_files * lf = &setting->local_files;        
    
    Kmer_handler kmer_handler(setting->kmer_length,  
            dup.min_count, dup.min_len, dup.threshold, dup.rmer_length, dup.is_use_set, lf);     
    
    Contig_handler contig_handler;       
    kmer_handler.get_contig_handler(&contig_handler);    
    
    if(kmer_handler.build_dict_from_kmer_file() != 0)
    {
        shc_log_error("Error reading kmers\n");        
    }
    kmer_handler.sort_kmer_descending_count(); 
    kmer_handler.find_contig(); 
    //kmer_handler.dump_kmers_to_file(lf->kmer_output_path);    
    //contig_handler.dump_all_contig(lf->contig_output_path);
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
