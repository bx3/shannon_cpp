#include "Contig_handler.h"

bool contig_handler_sort (contig_num_t i,contig_num_t j) { return (i<j); }

Contig_handler::Contig_handler(char * logname)
{
    num_contig = 0;
    shc_logname = logname;
    
    stage_count = 0;
    staged_num = 0x00;    
    delimitor.push_back(0);
}

Contig_handler::Contig_handler(char * logname, uint8_t rmer_length)
{
    num_contig = 0;
    shc_logname = logname;    
    stage_count = 0;
    staged_num = 0x00;
    delimitor.push_back(0);
}

void Contig_handler::dump_all_contig(std::string & filename)
{    
    shc_log_info(shc_logname, "Start Dumping all contigs into file\n");
    int reminder_bases_num = 0;
    std::ofstream outfile (filename.c_str());
    printf("numerb of contig is %u\n", num_contig);
    char four_base[4];
    for(int i=0; i<num_contig; i++)
    {       
        outfile << "Contig #" << i << " with mean count " << (int)mean_count[i] 
                << " : \n";
        // -1 since exclusive at the end
        //shc_log_info(shc_logname, "delimitor.at(i) : %d\n",delimitor.at(i));
        //shc_log_info(shc_logname, "delimitor.at(i+1) : %d\n",delimitor.at(i+1));
        for(int j=delimitor.at(i); j<delimitor.at(i+1); j++)  
        {
            if(contig_len_list[i]%FOUR_BASE==0 || j<delimitor.at(i+1)-1)
            {
                decode(contig_list.at(j), FOUR_BASE, four_base);
                for(int k=0; k<FOUR_BASE; k++)
                    outfile << four_base[k];
            }
            else
            {
                reminder_bases_num = contig_len_list[i]%FOUR_BASE;
                decode(contig_list.at(j), reminder_bases_num, four_base);
                for(int k=0; k<reminder_bases_num; k++)
                    outfile << four_base[k];
            }
        }
        outfile << std::endl;
    }
    shc_log_info(shc_logname, "Finish Dumping all contigs into file\n");
}

void Contig_handler::print_delimitor()
{ 
    for(std::vector<size_type>::iterator i = delimitor.begin(); i!=delimitor.end(); i++)
    {
        std::cout << *i << ' ';
    }
    std::cout << std::endl;
}

/**
 * Since we store all contig in a continuous vector, erasing elements helps 
 * to save storage, but we have to reallocate the following data's location
 * which is somehow inefficient, in the C++ standard, it takes linear time
 * 
 * Also notice to use C++ compiler, instead of the default one to use
 * shrink_to_fit function
 * @param i
 */
void Contig_handler::delete_contig(contig_num_t i)
{
    shc_log_info(shc_logname, "Start Deleting a single contig at index %u\n", i);
    contig_list.erase(contig_list.begin()+delimitor[i], 
                      contig_list.begin()+delimitor[i+1]); // since erase is [ )
    
    size_type num_base_deleted = delimitor[i+1] - delimitor[i];
    for(contig_num_t j=0; j<num_contig; j++)
    {
        if (j > i)
        {
            delimitor[j] = delimitor[j+1] - num_base_deleted;
        }
    }
    
    mean_count.erase(mean_count.begin()+i);
    num_contig--;    
    contig_list.shrink_to_fit();
    shc_log_info(shc_logname, "Finish Deleting a single contig at index %u\n", i);
}

/**
 * Providing a buffer for compressing bit
 * @param c
 */
void Contig_handler::push_back(uint8_t c)
{    
    staged_num = (staged_num<<2) | c;    
    shc_log_info(shc_logname, "staged %x with %c\n", staged_num, num_to_char[c]);
    if((++stage_count) >= FOUR_BASE)
    {
        shc_log_info(shc_logname, "Pushed %x\n", staged_num);
        contig_list.push_back(staged_num);
        stage_count = 0;
        staged_num = 0x00;        
    }        
}

/**
 * To flash the buffer in case no more input
 */
void Contig_handler::flash()
{
    contig_list.push_back(staged_num);
}

/**
 * Accepting this contig and update information
 * @param mean_c
 * @param contig_len
 */
void Contig_handler::declare_new_contig(kmer_count_t mean_c, size_type contig_len)
{    
    uint8_t * contig_start = &(contig_list.at(delimitor.back()));
    size_t byte_lenght = encode_in_place((char*)contig_start, contig_start, contig_len);
    //note the last is empty, and is reserved for delimitor to start
    contig_list.resize(delimitor.back()+byte_lenght); 
    num_contig++;
    delimitor.push_back(delimitor.back()+byte_lenght);
    mean_count.push_back(mean_c);
    contig_len_list.push_back(contig_len);   
    
    shc_log_info(shc_logname, "contig number increase to %u\n", num_contig-1);    
    shc_log_info(shc_logname, "delimitor adds an interval %u->%u \n", 
                               delimitor[delimitor.size()-2], delimitor[delimitor.size()-1]);            
    shc_log_info(shc_logname, "new contig has length %d\n", contig_len);      
    shc_log_info(shc_logname, "new contig has byte length %d\n", byte_lenght);      
}

/**
 * In case, contig is pushed in the list, but is not needed, resize the list
 * back to previous length.
 */
void Contig_handler::reject_new_contig()
{
    contig_list.resize(delimitor.back());
    shc_log_info(shc_logname, "reject contig, last delimitor remains %u\n", delimitor.back());    
}


/**
 * This function takes a list of contig index (contig number) and remove them
 * from the contig_list vector in a memory and algorithmic efficient way.
 * 
 * @param remove_list
 * @param remove_num
 */
void Contig_handler::delete_contig_list(contig_num_t * remove_list, contig_num_t remove_num)
{  
    shc_log_info(shc_logname, "Start Deleting a list of contigs\n");
    
    uint8_t * s = &contig_list[0];
    contig_num_t curr_index = remove_list[0];;          
             
    uint8_t * next = NULL;
    uint8_t * keep = NULL;
    uint8_t * mp = s + delimitor[curr_index];    
    
    size_type keep_index = 0;
    size_type last_to_keep_index = 0;
    
    size_type num_base_keep = 0;
    size_type num_base_remove = 0;
    size_type total_base_remove = 0;
    size_type keep_contig_num = 0;
    bool is_nothing_to_keep = false; 
    
    // in case we want to remove anything
    if(remove_num==num_contig)
    {
        std::vector<uint8_t>().swap(contig_list);
        num_contig = 0;
        std::vector<size_type>().swap(delimitor);
        std::vector<kmer_count_t>().swap(mean_count);
    }
        
    // going through the list and removing element, j is the removing index  
    // which points to the contig number to be removed
    for(contig_num_t j=0; j<remove_num; j++)
    {                
        curr_index = remove_list[j];        
                
        keep_index = get_next_to_keep(remove_list, remove_num, j);
        keep = s + delimitor[keep_index];
        
        if(keep_index == NOTHING_TO_KEEP )
        {
            is_nothing_to_keep = true;
            break;
        }
                
        contig_num_t next_delete_index = 
                get_next_to_delete(remove_list, remove_num, keep_index, j);                     

        next = s + delimitor[next_delete_index];
        keep_contig_num = next_delete_index - keep_index;
                                              
        // shifting the memory
        if( keep_index != last_to_keep_index)
        {     
            num_base_remove = delimitor[keep_index] - delimitor[curr_index]; 
            num_base_keep = delimitor[next_delete_index] - delimitor[keep_index];    
            
            total_base_remove += num_base_remove;
            std::cout << " curr_index " << curr_index << " is " << delimitor[curr_index] 
                      << " keep_index " << keep_index << " is " << delimitor[keep_index]
                      << " next_to_delete "<< next_delete_index    
                      << " num_base_remove "<< num_base_remove                       
                      << " num_base_keep "<< num_base_keep    
                      <<  std::endl;
                                    
            
                                
            memcpy(mp, keep, num_base_keep);
            last_to_keep_index = keep_index;                                         
            mp += num_base_keep;         
        }        
    }  
    
    std::cout << "total_base_remove " << total_base_remove<< std::endl;
    
    if(is_nothing_to_keep)
    {        
        
        contig_list.resize(delimitor[curr_index] - total_base_remove);
        contig_list.shrink_to_fit();           
    }
    else
    {                
        contig_list.resize(contig_list.size()-total_base_remove);
        contig_list.shrink_to_fit();   
    }
        
    update_delimitor_and_mean_count(remove_list,  remove_num);    
    num_contig -= remove_num;        
    delimitor.resize(num_contig+1);
    delimitor.shrink_to_fit(); 
    shc_log_info(shc_logname, "Finish Deleting a list of contigs\n");
}

/**
 * This function uses the current remove contig index to locate the next 
 * available index of contig to keep. remove list is sorted in ascending order
 * @param remove_list
 * @param remove_num
 * @param curr
 * @return 
 */

contig_num_t Contig_handler::get_next_to_keep(
                        contig_num_t * remove_list, contig_num_t remove_num, 
                        contig_num_t remove_list_index)
{  
    contig_num_t remove_index;
    
    for(contig_num_t i=remove_list_index; i<remove_num; i++)
    {   
        remove_index = remove_list[i];
        
        //last item in remove_list
        if(i==remove_num-1)
        {
            if (remove_index == num_contig-1)
                return NOTHING_TO_KEEP;
            return remove_index+1;
        }
        
        //check if next entry available
        if(remove_list[i+1] != remove_index+1)
        {
            return remove_index+1;
        }                          
    }        
}
            
contig_num_t Contig_handler::get_next_to_delete(
        contig_num_t * remove_list,  contig_num_t remove_num, 
        contig_num_t keep_index, contig_num_t remove_index)
{ 
    if(remove_num==1)
        return num_contig;
    
    contig_num_t last_index = remove_num-1;
    contig_num_t last = remove_list[last_index];
    contig_num_t i = 0;
    for(i=remove_num-1; i>0; i--)
    {        
        last = remove_list[i];
        std::cout << "last " <<last << std::endl;
        if(last-1 != remove_list[i-1])
        {
            last_index = i;
            break;
        }
    }
    if(i==0)
    {
        return num_contig;
    }
            
    std::cout << "last_index " <<last_index << std::endl;
    std::cout << "remove_index " <<remove_index << std::endl;
    
    contig_num_t counter = 1;
    contig_num_t next_remove = remove_list[remove_index+counter];
    for(contig_num_t i=remove_list[remove_index]+1; i<num_contig; i++)
    {
        if(remove_index == remove_num-1)
        {
            return num_contig;
        }
        
        if(remove_index == last_index)
        {
            return num_contig;
        }
                
        if(next_remove == i)
        {
            counter++;
            next_remove = remove_list[counter+remove_index];
        }
        else
        {
            std::cout << "next_remove " <<next_remove << std::endl;
            return next_remove;
        }
    }
}


void Contig_handler::print_contig_length()
{
    for(size_type i=0; i< num_contig; i++)
        std::cout << delimitor[i+1] - delimitor[i] << " ";
    std::cout << std::endl;    
}

void Contig_handler::print_delimitor(std::vector<size_type> & d)
{
    for(std::vector<size_type>::iterator i = d.begin(); i!=d.end(); i++)
    {
        std::cout << *i << ' ';
    }
    std::cout << std::endl;
}

void Contig_handler::print_contig(size_type contig_start, size_type base_len)
{
    shc_log_info(shc_logname, "contig length %d\n", base_len);    
    uint8_t re = contig_list[contig_start+base_len];
    contig_list[contig_start+base_len] = '\0';
    uint8_t * print_start = &(contig_list.at(contig_start));
    shc_log_info(shc_logname, "contig!!: %s\n", (char*)print_start);    
    contig_list[contig_start+base_len] = re;
}

void Contig_handler::update_delimitor_and_mean_count(contig_num_t * remove_list, contig_num_t remove_num)
{
    //construct the list
    std::vector<size_type> length_list;    
    std::vector<size_type> new_length_list;  
    std::vector<size_type> new_delimitor;  
    std::vector<kmer_count_t> new_mean_count;
    
    for(size_type i=0; i<num_contig; i++)
    {
        length_list.push_back(delimitor[i+1]-delimitor[i]);        
    }
    
    // remove the contig in length array
    contig_num_t remove_counter = 0;
    contig_num_t coming_remove_contig = remove_list[remove_counter];    
    for(contig_num_t i=0; i<num_contig; i++)
    {
        //std::cout << "coming_remove_contig is " << coming_remove_contig<< std::endl;
        if (i!=coming_remove_contig)
        {
            new_length_list.push_back(length_list[i]);    
            new_mean_count.push_back(mean_count[i]);
        }
        else
        { 
            ++remove_counter;
            if(remove_counter<remove_num)            
                coming_remove_contig = remove_list[remove_counter];                 
            else
                coming_remove_contig = num_contig;            
        }                
    }
    
    new_delimitor.push_back(0);
    for(std::vector<size_type>::iterator i=new_length_list.begin();
            i<new_length_list.end(); i++)
    {
        new_delimitor.push_back(new_delimitor.back()+ *i );
    }
    delimitor.swap(new_delimitor);
    delimitor.resize(num_contig+1-remove_num);
    delimitor.shrink_to_fit();
    
    mean_count.swap(new_mean_count);
    mean_count.resize(num_contig-remove_num);
    mean_count.shrink_to_fit();
    
    //free memory
    std::vector<size_type>().swap(new_delimitor);
    std::vector<size_type>().swap(new_length_list);
    std::vector<size_type>().swap(length_list);
    std::vector<kmer_count_t>().swap(new_mean_count);
}

void Contig_handler::get_logname(char * logname)
{
    shc_logname = logname;
}

/*
void Contig_handler::filter_contig(kmer_count_t min_count, size_type min_len,
                                   double threshold)
{    
    shc_log_info(shc_logname, "Start Filtering contigs\n");
    shc_log_info(shc_logname, "Start length count checks\n");
    std::vector<contig_num_t> remove_list;
    for(contig_num_t i=0; i<num_contig; i++)
    {
        size_type len = delimitor[i+1] - delimitor[i];
        if(len<min_len || mean_count[i]<min_count)
        {
            shc_log_info(shc_logname, "Remove %u \n", i);
            remove_list.push_back(i);
        }
    }        
                
    if(!remove_list.empty())
    {
        get_vector_string(remove_list, vec_str);
        shc_log_info(shc_logname, "Deleting vec %s\n", vec_str.str().c_str());
        delete_contig_list( &(remove_list.at(0)), remove_list.size());            
    }
    else    
    {
        shc_log_info(shc_logname, "length count checks and Nothing to delete\n");    
    }
    
    remove_redundance(threshold);
    
    std::vector<contig_num_t>().swap(remove_list);
    shc_log_info(shc_logname, "Finish Filtering contigs\n");
}

void Contig_handler::remove_redundance(double threshold)
{
    shc_log_info(shc_logname, "Start Removing redundance with total %u contigs\n", num_contig);
    size_type contig_begin = 0;
    size_type contig_end = 0;
    
    uint8_t * rmer_begin = NULL;    
    size_type contig_size = 0;
    char temp[30];
    temp[rmer_len] = '\0';
    
    double i_thresh, j_thresh;
    
    uint64_t byte = 0;
    // for each contig, i which contig
    for(contig_num_t i=0; i<num_contig; i++)
    {
        contig_begin = delimitor[i];
        contig_end = delimitor[i+1];
        contig_size = contig_end - contig_begin;
        // get its k-1 kmer
        for(size_type j=0; j<contig_size - rmer_len ; j++)
        {
            rmer_begin = &(contig_list.at(contig_begin + j));            
            encode_kmer((char*)rmer_begin, &byte, rmer_len);
            //exclude the case when r are made of the same character, and the
            //corresponding repeat in the vector
            if (rmer_contig_map[byte].empty())
            {
                rmer_contig_map[byte].push_back(i);            
            }
            else if(!rmer_contig_map[byte].empty() && rmer_contig_map[byte].back()!= i)
            {   
                 // debug
                //memcpy(temp, rmer_begin, rmer_len);
                //printf("contig %u:%u: %s\n", i, j, temp);                                            
                rmer_contig_map[byte].push_back(i);
            }                 
        }
    }
    
    std::vector<contig_num_t> remove_indicator(num_contig, 0);
    std::vector<contig_num_t> remove_list;
    size_type i_begin, i_end, i_size, j_size;   
    uint64_t counter = 0;
    bool is_find = false;
    for(contig_num_t ci=num_contig; ci!=0; ci--)
    {
        contig_num_t i = ci - 1;
        shc_log_info(shc_logname, "Redundance examing %u\n", i);
        if(remove_indicator.at(i) == 0)
        {
            i_begin = delimitor[i];
            i_end = delimitor[i+1];
            i_size = i_end - i_begin;            
            //compare with other contig
            for(contig_num_t cj=ci-1; cj!=0; cj--)
            {
                contig_num_t j = cj - 1;                
                j_size = delimitor[j+1] - delimitor[j];
                counter = 0;
                // find its little R-mer             
                for(size_type j=0; j<contig_size - rmer_len ; j++)
                {
                    rmer_begin = &(contig_list.at(contig_begin + j));            
                    encode_kmer((char*)rmer_begin, &byte, rmer_len);

                    is_find = std::binary_search(rmer_contig_map[byte].begin(), 
                                                 rmer_contig_map[byte].end(), 
                                                 j);
                    if(is_find)                    
                        counter++;
                }
                
                i_thresh = (double)(counter*1.0/(i_size));
                j_thresh = (double)(counter*1.0/(j_size));
                
                shc_log_info(shc_logname, "comparing %u -> %f with %u-> %f, count %u %u\n", 
                                                     i, i_thresh,   j,  j_thresh, 
                                                     i_size, j_size);
                if(i_size > j_size)
                {
                    if(j_thresh > threshold)
                    {
                        remove_indicator[j] = 1;  
                        remove_list.push_back(j);
                        
                        get_vector_string(remove_list, vec_str);
                        shc_log_info(shc_logname, "Deleting vec %s\n", vec_str.str().c_str());                                        
                        
                    }
                    
                }
                else
                {
                    if(i_thresh > threshold)
                    {
                        remove_indicator[i] = 1;                                        
                        remove_list.push_back(i);
                        
                        get_vector_string(remove_list, vec_str);
                        shc_log_info(shc_logname, "Deleting vec %s\n", vec_str.str().c_str());                                        
                        
                        break;
                    }                    
                }                           
            }
        }
    }
    if(!remove_list.empty())
    {
        get_vector_string(remove_list, vec_str);
        shc_log_info(shc_logname, "Redundance checks to delete %s\n", vec_str.str().c_str());
            
        std::sort(remove_list.begin(), remove_list.end(), contig_handler_sort);
        
        get_vector_string(remove_list, vec_str);
        shc_log_info(shc_logname, "Redundance checks to delete sorted %s\n", vec_str.str().c_str());
        //free memory
        delete_contig_list( &(remove_indicator.at(0)), remove_indicator.size());    
        Rmer_counter_map().swap(rmer_contig_map);
    }
    else
    {
        shc_log_info(shc_logname, "Redundance checks and Nothing to delete\n");
    }
    std::vector<contig_num_t>().swap(remove_indicator);
    
    shc_log_info(shc_logname, "Finishing Removing redundance\n");
}
 */