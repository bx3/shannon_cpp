#ifndef MULTIBRIDGE_HANDLER_H
#define MULTIBRIDGE_HANDLER_H

#include <pthread.h>
#include "Sequence_graph_handler.h"
#include "Kmer_handler.h"
#include "shc_type.h"
#include "json_parser.h"
#include "local_file_structure.h"
#include <vector>
#include <list>

class Multibridge_handler {
    friend class Contig_handler;
public:
    Multibridge_handler(Shannon_C_setting & setting_, uint8_t kmer_length_,
                        bool is_compressed_) : 
                        setting(setting_), kmer_length(kmer_length_), 
                        is_compressed(is_compressed_) {}
    Multibridge_handler(Kmer_handler * kh_, Shannon_C_setting & setting_, 
                        uint8_t kmer_length_, bool is_compressed_) : 
                         kh(kh_), setting(setting_), kmer_length(kmer_length_), 
                        is_compressed(is_compressed_) {}
    
    void count_num_component_with_file(std::string & read_dir, std::string & kmer_dir);
    
    void process_components(int num_thread);
    
    
        
private:
    std::list<Sequence_graph_handler *> comp_graph_list;     
    int num_comp;
    uint8_t kmer_length;
    bool is_compressed;
    Kmer_handler * kh;   
    Shannon_C_setting setting;
    
    Block_timer timer;
};

#endif

