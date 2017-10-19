#include "Multibridge_handler.h"

void Multibridge_handler::
count_num_component_with_file(std::string & read_path, std::string & kmer_path)
{   
    int read_counter=0, kmer_counter=0;          
        
    for(boost::filesystem::directory_iterator it(read_path);
        it != boost::filesystem::directory_iterator(); it++)
    {
        read_counter++;
    }
    
    for(boost::filesystem::directory_iterator it(kmer_path);
        it != boost::filesystem::directory_iterator(); it++)
    {
        kmer_counter++;
    }
    assert(read_counter == kmer_counter);    
    num_comp = read_counter;
    std::cout << "number of component is " << num_comp << std::endl;
}

void Multibridge_handler::process_components(int num_thread)
{
    start_timer(&timer);
    void * tret;
    int max_hop = 100000;
    std::list<pthread_t> threads;
    if(num_comp < num_thread)
    {   
        int i = 0;
        while(i < num_comp)
        {
            Sequence_graph_handler * seq_graph_handler = new 
                    Sequence_graph_handler(kmer_length, setting, i, is_compressed, max_hop);
            comp_graph_list.push_back(seq_graph_handler);

            pthread_t thread;
            
            int status = pthread_create(&thread, NULL, 
                seq_graph_handler->run_seq_graph_handler, (void*)seq_graph_handler);
            if(status != 0)
            {
                std::cerr << "thread creation fails" << std::endl;
                exit(1);
            }
            else
            {
                threads.push_back(thread);
                std::cout << "created thread " << i << std::endl;
                i++;
            }
        }
    }
    else
    {
        std::cerr << "case not handled" << std::endl;
        exit(1);
    }
    
    for(std::list<pthread_t>::iterator it=threads.begin();
                                       it!=threads.end(); ++it)
    {
        //std::cout << "before join" << std::endl;        
        pthread_join(*it, &tret);
        //std::cout << "after join" << std::endl;        
    }
    
    std::cout << "multibridge finish, ";
    stop_timer(&timer);
}