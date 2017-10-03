#include "Contig_graph_handler.h"

size_t contig_curr_proc = 0;
size_t contig_prev_proc = 0;

struct Block_timer cgh_timer;

void Contig_graph_handler::group_components()
{          
    shc_log_info(shc_logname, "Start constructing contig graph\n");      
    Kmer_counter_map_iterator it;
    char kmer_array[33];
    kmer_array[kh->kmer_length] = '\0';
    uint64_t byte = 0;    
            
    //put every contig into the set
    for(contig_num_t i=0; i<ch->num_contig; i++)
        explorable_contig_set.insert(i);   
    
    start_timer(&cgh_timer);
    while(!explorable_contig_set.empty())
    {
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "\t\t\t\t\t#processing component %u\n",curr_component_num);
#endif
#ifdef SHOW_PROGRESS
        if((contig_curr_proc-contig_prev_proc) > CONTIG_PROGRESS_STEP)
        {
            int percentage = (100 * contig_curr_proc)/(ch->num_contig);
            std::cout << "[" << percentage << "%] " <<  contig_curr_proc 
                 << " contig out of " <<   ch->num_contig
                 << " is processed" << std::endl;
            contig_prev_proc = contig_curr_proc;
        }
#endif        
        
        contig_num_t root_contig = *explorable_contig_set.begin();
        contig_stack.push(root_contig);
        explorable_contig_set.erase(explorable_contig_set.begin());
        
        // add that node, so each graph has at least one node
        if(contig_vertex_map.find(root_contig) == contig_vertex_map.end())
        {
            vd_t vd = boost::add_vertex(graph);
            contig_vertex_map[root_contig] = vd;            
            graph[vd] = bundled_contig_index(root_contig);
            //shc_log_info(shc_logname, "Root contig %u with vd %u\n", root_contig, vd);
        }
        
        //explore that contig and find new contig
        while(!contig_stack.empty())
        {
            contig_num_t i = contig_stack.top();
            contig_stack.pop();
            
            Contig_handler::Contig_base_info contig_info = 
                    ch->get_contig(i);
            
            for(Contig_handler::size_type j=0;
                j<ch->contig_len_list.at(i)-kh->kmer_length+1;   j++)
            {                
                memcpy(kmer_array, contig_info.base_start+j, kh->kmer_length);
                encode_kmer(contig_info.base_start+j, &byte, kh->kmer_length);
                
                //shc_log_info(shc_logname, "%s\n", kmer_array);
                
                it = kh->kmer_counter.find(byte);                                
                if(it == kh->kmer_counter.end())
                {                    
                    shc_log_error("contain impossible kmer %s\n", kmer_array);
                    exit(1);
                }                
                get_connected_contig_index(it);
            }            
        }
        contig_curr_proc = ch->num_contig - explorable_contig_set.size();
        //up to here we have a component and its graph
        break_and_keep_component();
        graph.clear();
        contig_vertex_map.clear();
    }
    
#ifdef SHOW_PROGRESS
    std::cout << "Finish finding component, ";
    stop_timer(&cgh_timer);
#endif
    
    shc_log_info(shc_logname, "Finish constructing contig graph\n");
}

void Contig_graph_handler::break_and_keep_component()
{   
    typedef std::pair<vi_t, vi_t> vip_t;
    idx_t num_vertices = boost::num_vertices(graph);
    //shc_log_info(shc_logname, "num vertices %u, partition_size %u\n",num_vertices , metis_setup.partition_size);
    if( num_vertices <= metis_setup.partition_size)
    {
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "no metis\n");
#endif
        if(!is_set_collect_comp_num)
        {
            collect_comp_num = curr_component_num;
            curr_component_num++;    
            is_set_collect_comp_num = true;
        }
        
        vip_t vip = vertices(graph);
        //update contig to component info
        for(vi_t it=vip.first; it!=vip.second; it++)
        {
            component_array[graph[*it].contig_index] = collect_comp_num;  
            //shc_log_info(shc_logname, "component %u has contig %u\n", curr_component_num, graph[*it].contig_index);
        }          
    }
    else        //call metis
    {
#ifdef LOG_METIS
        shc_log_info(shc_logname, "Starting metis \n");
#endif
        create_metis_array_input();
        metis_setup.num_vertices = num_vertices;
        metis_setup.num_partition = std::min(static_cast<idx_t>(100), 
                num_vertices/metis_setup.partition_size + 1);
        //shc_log_info(shc_logname, "metis to asked to break %d partition\n", metis_setup.num_partition);
        idx_t ufactor = static_cast<idx_t>(METIS_IMBALANCE_PRECISION 
                                                   * (metis_setup.overload-1));  
        metis_setup.options[METIS_OPTION_UFACTOR] = ufactor; 
        metis_setup.options[METIS_OPTION_NUMBERING] = 0;
        
        std::vector<idx_t> part(metis_setup.num_vertices);
        
        run_single_pass_metis(component_array, &part);
        if(metis_setup.is_multiple_partition)
        {
            add_weight_to_cut_and_update_graph(&part);            
            run_single_pass_metis(component_array_aux, &part);
        }
        std::vector<idx_t>().swap(part);    
    }
#ifdef LOG_CONTIG_GRAPH
    shc_log_info(shc_logname, "end\n");
#endif    
    
}

void Contig_graph_handler::run_single_pass_metis(std::vector<comp_num_t> & array, std::vector<idx_t> * part)
{                             
    int ret = METIS_PartGraphRecursive(
          &metis_setup.num_vertices, &metis_setup.ncon, 
          &metis_input.xadj.at(0), &metis_input.adjncy.at(0),
          NULL, NULL, &metis_input.adjwgt.at(0), &metis_setup.num_partition, 
          NULL, NULL, NULL, &metis_setup.objval, &part->at(0));
    
    if(ret == METIS_ERROR_INPUT)
        shc_log_error("METIS_ERROR_INPUT\n");
    if(ret == METIS_ERROR_MEMORY)
        shc_log_error("METIS_ERROR_MEMORY\n");
    if(ret == METIS_ERROR)
        shc_log_error("METIS_OTHER_ERROR\n");
        
#ifdef LOG_METIS    
    shc_log_info(shc_logname, "metis breaks to %d partition\n", metis_setup.num_partition);
#endif        
        
    for(idx_t i=0; i<metis_setup.num_vertices; i++)
    {
        array[graph[static_cast<vd_t>(i)].contig_index] = 
                                             curr_component_num + part->at(i);
    }
    
    curr_component_num += metis_setup.num_partition;    
}

void Contig_graph_handler::add_weight_to_cut_and_update_graph(std::vector<idx_t> * part)
{
    typename boost::graph_traits < graph_t >::adjacency_iterator vi, vi_end; 
    std::pair<ed_t, bool> aer;
    std::pair<vi_t, vi_t> vip = vertices(graph);
    //walk through adj list, and update if 
    for(vi_t it=vip.first; it!=vip.second; it++)
    {        
        vd_t vd_source = *it; 
        idx_t source_comp = part->at(static_cast<idx_t>(vd_source));
        for (boost::tie(vi, vi_end)=boost::adjacent_vertices(vd_source, graph); 
                vi != vi_end; ++vi)
        {
            idx_t adj_comp = part->at(static_cast<idx_t>(*vi));
            if(adj_comp != source_comp)
            {
                aer = boost::edge(vd_source, *vi , graph);  
                graph[aer.first].weight *= metis_setup.penalty;
                //shc_log_info(shc_logname, "edge between %d and %d is increased to %d\n",vd_source, *vi, graph[aer.first].weight);
            }
        }                       
    }
    create_metis_array_input();
}

/*
if(aer.second)
    {                
        graph[aer.first].weight++;
    }
*/

void Contig_graph_handler::get_connected_contig_index(
                                                Kmer_counter_map_iterator & it)
{        
    contig_num_t contigA_index, contigB_index;
    uint64_t next_byte = 0;    
    contigA_index = it->second.contig;    
    // for that kmer, check what extension available to it
    for(uint8_t i=1; i<= BIT_PER_BYTE; i++)
    {
        if(kh->is_info_ith_bit_set(it->second.info, i))
        {
            //shc_log_info(shc_logname, "%lld : %u bit is set \n",it->first, i);
            if(i<=PREFIX_OFFSET)           
                next_byte = prepend_byte(&(it->first), i-1, kh->kmer_length);
            else            
                next_byte = append_byte(&(it->first), i-1-PREFIX_OFFSET, kh->kmer_length);                     
            //in case the matched kmer is deleted since that does not belong 
            //to any contig or that contig is deleted
            if (kh->kmer_counter.find(next_byte) != kh->kmer_counter.end())
            {
                contigB_index = kh->kmer_counter[next_byte].contig;
                if (contigA_index != contigB_index)
                {    
                    //shc_log_info(shc_logname, "Find new contig %u, before increment edge with me %u\n", contigB_index, contigA_index);
                    increment_edge_weight(contigA_index, contigB_index);                    
                    //shc_log_info(shc_logname, "AFter increment");
                    if (explorable_contig_set.find(contigB_index) != 
                        explorable_contig_set.end())
                    {
                        //if contigB_index is found
                        contig_stack.push(contigB_index);
                        explorable_contig_set.erase(contigB_index);
                        //shc_log_info(shc_logname, "%u is pushed to search stack\n", contigB_index);
                    }
                }
            }
        }
    }
}

/**
 * Add an edge between two nodes, if there isn't a node, add an edge, otherwise
 * increment by 1, since edge is number of common k-1 mer
 * @param i
 * @param j
 */
void Contig_graph_handler::increment_edge_weight(contig_num_t i, contig_num_t j)
{    
    //add nodes if needed    
    vd_t vd;
    bool is_wihtout = contig_vertex_map.find(i) == contig_vertex_map.end();    
    if(is_wihtout)
    {
        //shc_log_info(shc_logname, "without %u\n", i);
        vd = boost::add_vertex(graph);
        contig_vertex_map[i] = vd;   
        graph[vd] = bundled_contig_index(i);
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "Added contig %u with vd %u\n", i, vd);
#endif
    }
    is_wihtout = contig_vertex_map.find(j) == contig_vertex_map.end();
    
    if(is_wihtout)
    {
        //shc_log_info(shc_logname, "without %u\n", j);
        vd = boost::add_vertex(graph);
        contig_vertex_map[j] = vd;        
        graph[vd] = bundled_contig_index(j);
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "Added contig %u with vd %u\n", j, vd);
#endif
    }
        
    std::pair<ed_t, bool> aer;
        
    aer = boost::edge(contig_vertex_map[i], contig_vertex_map[j], graph);    
    if(aer.second)
    {                
        graph[aer.first].weight++;
        
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "edge between %u %u : weight increased to %u\n",
                i, j ,graph[aer.first].weight);        
#endif
    }
    else //add edge
    {
        aer = boost::add_edge(contig_vertex_map[i], contig_vertex_map[j], graph);
        assert(aer.second);
        graph[aer.first] = bundled_weight(1);                
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "Edge between %u %u\n", i, j);
#endif
    }            
}

void Contig_graph_handler::create_metis_array_input()
{
    shc_log_info(shc_logname, "Start creating metis array input for component\n",
                                                curr_component_num);    
    typename boost::graph_traits < graph_t >::adjacency_iterator vi, vi_end;    
    vd_t vd;    
    std::pair<ed_t, bool> aer;
    contig_edge_weight_t edge_weight;
    
    auto vip = vertices(graph);
    
    for(vi_t it=vip.first; it!=vip.second; it++)
    {        
        vd_t vd_source = *it;                
        for (boost::tie(vi, vi_end)=boost::adjacent_vertices(vd_source, graph); 
                vi != vi_end; ++vi)
        {
            vd = *vi;
            aer = boost::edge(vd_source, vd, graph);
            // since for AAAT, AATG we double counts
            edge_weight = graph[aer.first].weight/2; 
            metis_input.adjncy.push_back(vd);
            metis_input.adjwgt.push_back(edge_weight);            
        }  
        metis_input.xadj.push_back(metis_input.adjncy.size());                
    }    
    //log_metis_input_data();
    shc_log_info(shc_logname, "Finish creating metis array input for component\n");
}

void Contig_graph_handler::create_metis_file_format_from_graph(std::string & filename)
{
    shc_log_info(shc_logname, "Start making METIS format file\n");    
    typename boost::graph_traits < graph_t >::adjacency_iterator vi, vi_end;    
    vd_t vd;
    std::ofstream outfile (filename);
    std::pair<ed_t, bool> aer;
    contig_edge_weight_t edge_weight;
    //first line in metis file
    outfile << boost::num_vertices(graph) << " " << boost::num_edges(graph) 
            << " " << ONLY_EDGE_WEIGHT << std::endl;       
            
    auto vip = vertices(graph);
    
    for(vi_t it=vip.first; it!=vip.second; it++)
    {        
        vd_t vd_source = *it;        
        //uncomment if want to see the graph structure in file clear
        //outfile << vd_source << ":\t" ;
        for (boost::tie(vi, vi_end)=boost::adjacent_vertices(vd_source, graph); 
                vi != vi_end; ++vi)
        {
            vd = *vi;
            aer = boost::edge(vd_source, vd, graph);
            // since for AAAT, AATG we double counts
            edge_weight = graph[aer.first].weight/2; 
            metis_input.adjncy.push_back(vd);
            metis_input.adjwgt.push_back(edge_weight);
            //outfile << vd+METIS_NODE_OFFSET << " " << edge_weight << " ";
            outfile << vd << " " << edge_weight << " ";    //without offset               
        }  
        metis_input.xadj.push_back(metis_input.adjncy.size());
        
        outfile << std::endl;
    }
    
    //log_metis_input_data();
    outfile.close();
    shc_log_info(shc_logname, "Finish making METIS format file\n");
} 


void Contig_graph_handler::log_metis_input_data()
{
    shc_log_info(shc_logname, "Printing metis input data\n");
    
    info_log_str_without_new_line(shc_logname,"\n xadj:   ");
    for(auto iter=metis_input.xadj.begin(); iter!=metis_input.xadj.end(); ++iter)    
        info_log_num_without_new_line(shc_logname, *iter);        
        
    info_log_str_without_new_line(shc_logname,"\n adjncy: ");
    for(auto iter=metis_input.adjncy.begin(); iter!=metis_input.adjncy.end(); ++iter)    
        info_log_num_without_new_line(shc_logname, *iter);        
    
    info_log_str_without_new_line(shc_logname,"\n adjwgt: ");
    for(auto iter=metis_input.adjwgt.begin(); iter!=metis_input.adjwgt.end(); ++iter)    
        info_log_num_without_new_line(shc_logname, *iter);     
    info_log_str_without_new_line(shc_logname,"\n\n");
}

void Contig_graph_handler::dump_component_array(std::string & filename)
{    
    shc_log_info(shc_logname, "Start dumping: %u component \n", curr_component_num);
    std::ofstream outfile(filename);    
    outfile << curr_component_num << " components" << std::endl;
    outfile << "contig"  << "\t"  << "components" << std::endl;
    //int i=0;
    for(int i=0; i<component_array.size(); i++)
    {
        if(!metis_setup.is_multiple_partition)
        {
            outfile << component_array[i] << std::endl;
        }
        else
        {
            outfile << component_array[i] << "\t" << component_array_aux[i] << std::endl;
        }
        //shc_log_info(shc_logname, "%d   %u\n", i++, *it);
    }
    outfile.close();
    shc_log_info(shc_logname, "Finish log_component_array\n");
}

void Contig_graph_handler::assign_paired_read_to_components(int num_test)
{  
    
    shc_log_info(shc_logname, "Start assigning read to component\n");
    //creating dir
    typedef boost::filesystem::path path_t;        
    std::string dir_path = lf->output_path_str + "/components_reads";  
    path_t comp_read_path(dir_path);
    if(boost::filesystem::exists(comp_read_path))
        replace_directory(comp_read_path);
    else
        add_directory(comp_read_path);

    //create and open files
    std::vector<std::shared_ptr<std::ofstream> > files;        
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        std::string file_path = dir_path + 
                        (std::string("/comp") + std::to_string(i));
        std::shared_ptr<std::ofstream> file(new std::ofstream);
        file->open(file_path.c_str());        
        files.push_back(file);
    }
    shc_log_info(shc_logname, "Created contig files\n"); 

    std::ifstream file1_reader(lf->input_read_path.c_str());  
    std::ifstream file2_reader(lf->input_read_path_2.c_str());  
    shc_log_info(shc_logname, "successfully open read files\n");
    std::string line1, header1, sequence1; 
    std::string line2, header2, sequence2;     
    int interval1 = 1, interval2 = 1;
    uint64_t byte1 = 0, byte2 = 0;    
    
    std::map<comp_num_t, int> comp_count;
    std::map<comp_num_t, int> comp_count_aux;
        
    while (file1_reader.good() && file2_reader.good() ) 
    {
        std::getline(file1_reader, line1);
        std::getline(file2_reader, line2);
        //std::cout << "content " << line1 << " " << line2 << std::endl;
        if(line1[0] == '>' && line2[0] == '>')
        {                      
            if(!sequence1.empty() && !sequence2.empty())
            {
                //process                
                interval1 = (sequence1.size()-kh->kmer_length+1)/num_test;
                interval2 = (sequence2.size()-kh->kmer_length+1)/num_test;
                for(int i=0; i<num_test; i++)
                {
                    //shc_log_info(shc_logname, "Before encoding readsize %u at %d\n",sequence.size(), i*interval);                                        
                    encode_kmer(&sequence1.at(0)+i*interval1, &byte1, kh->kmer_length);
                    encode_kmer(&sequence2.at(0)+i*interval2, &byte2, kh->kmer_length);
                    Kmer_counter_map_iterator it1 = kh->kmer_counter.find(byte1);
                    Kmer_counter_map_iterator it2 = kh->kmer_counter.find(byte2);
                    
                    if(it1!= kh->kmer_counter.end())
                    {                            
                        
                        contig_num_t contig_num_for_kmer = (kh->kmer_counter[byte1]).contig;                        
                        comp_count[component_array[contig_num_for_kmer]]++;            
                        
                        if(metis_setup.is_multiple_partition && 
                           component_array_aux[contig_num_for_kmer] < curr_component_num
                           )                            
                            comp_count_aux[component_array_aux[contig_num_for_kmer]]++;
                        
                        
                    }                    
                    if(it2!= kh->kmer_counter.end())
                    {                        
                        contig_num_t contig_num_for_kmer = (kh->kmer_counter[byte1]).contig;                        
                        comp_count[component_array[contig_num_for_kmer]]++;            
                        
                        if(metis_setup.is_multiple_partition && 
                           component_array_aux[contig_num_for_kmer] < curr_component_num
                           )                            
                            comp_count_aux[component_array_aux[contig_num_for_kmer]]++;
                    }
                }     
                // check if this sequence
                if(!comp_count.empty())                
                {
                    auto max_iter = std::max_element(comp_count.begin(), comp_count.end(), 
                            [](const std::pair<comp_num_t, int> & p1, 
                               const std::pair<comp_num_t, int> & p2 )
                            {
                                return p1.second < p2.second;
                            }
                    );

                    comp_num_t best_comp =  max_iter->first; 
                    //shc_log_info(shc_logname, "Before add to file %s \n",sequence.c_str());
                    *(files.at(best_comp)) << header1 << std::endl;
                    *(files.at(best_comp)) << sequence1 << std::endl;
                    *(files.at(best_comp)) << header2 << std::endl;
                    *(files.at(best_comp)) << sequence2 << std::endl;
                    //shc_log_info(shc_logname, "after add to file\n");                                
                    comp_count.clear();                                
                    sequence1.clear();
                    sequence2.clear();
                }
                
                if(metis_setup.is_multiple_partition && !comp_count_aux.empty())                
                {
                    auto max_iter = std::max_element(comp_count_aux.begin(), comp_count_aux.end(), 
                            [](const std::pair<comp_num_t, int> & p1, 
                               const std::pair<comp_num_t, int> & p2 )
                            {
                                return p1.second < p2.second;
                            }
                    );

                    comp_num_t best_comp =  max_iter->first; 
                    
                    //shc_log_info(shc_logname, "Before add to file %s \n",sequence.c_str());
                    *(files.at(best_comp)) << header1 << std::endl;
                    *(files.at(best_comp)) << sequence1 << std::endl;
                    *(files.at(best_comp)) << header2 << std::endl;
                    *(files.at(best_comp)) << sequence2 << std::endl;
                    //shc_log_info(shc_logname, "after add to file\n");                                
                    comp_count_aux.clear();                                                    
                }
            }
            header1 = line1;
            header2 = line2; 
            //std::cout << header1 << " " << header2 << std::endl;
        }
        else //line is sequence
        {     
            sequence1 += line1;
            sequence2 += line2;
        }                          
    }
    
    
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        (*files.at(i)).close();                
    }
    
    std::vector<std::shared_ptr<std::ofstream> >().swap(files);
    shc_log_info(shc_logname, "Finish assigning read to component\n");
    
}


void Contig_graph_handler::assign_reads_to_components(int num_test )
{
    //Fasta_reader fasta_reader("SE_read.fasta");        
    shc_log_info(shc_logname, "Start assigning read to component\n");
    //creating dir
    typedef boost::filesystem::path path_t;        
    
    std::string dir_path = lf->output_path_str + "/components_reads";  
    path_t comp_read_path(dir_path);
    if(boost::filesystem::exists(comp_read_path))
        replace_directory(comp_read_path);
    else
        add_directory(comp_read_path);
    shc_log_info(shc_logname, "Created dir %s\n", dir_path.c_str());
    std::cout << "created dir " << std::endl;
            
    //create files    
    std::vector<std::shared_ptr<std::ofstream> > files;
    
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        std::string file_path = dir_path + 
                        (std::string("/comp") + std::to_string(i));
        std::shared_ptr<std::ofstream> file(new std::ofstream);
        file->open(file_path.c_str());
        
        files.push_back(file);
    }
    shc_log_info(shc_logname, "Created contig files\n");
    
    //read and process files
    std::ifstream file_reader(lf->input_read_path);  
    shc_log_info(shc_logname, "successfully open read files\n");
    std::string line;
    std::string header;
    std::string sequence;    
    int interval = 1;
    uint64_t byte = 0;    
    
    std::map<comp_num_t, int> comp_count;
    std::map<comp_num_t, int> comp_count_aux;
        
    while (file_reader.good() ) 
    {
        std::getline(file_reader, line);
        if(line[0] == '>')
        {                     
            if(!sequence.empty())
            {
                //process                
                interval = (sequence.size()-kh->kmer_length+1)/num_test;
                for(int i=0; i<num_test; i++)
                {
                    //shc_log_info(shc_logname, "Before encoding readsize %u at %d\n",sequence.size(), i*interval);                                        
                    encode_kmer(&sequence.at(0)+i*interval, &byte, kh->kmer_length);
                    Kmer_counter_map_iterator it = kh->kmer_counter.find(byte);
                    if(it!= kh->kmer_counter.end())
                    {
                        contig_num_t contig_num_for_kmer = (kh->kmer_counter[byte]).contig;
                        //shc_log_info(shc_logname, "contig_num_for_kmer %u\n", contig_num_for_kmer );
                        //shc_log_info(shc_logname, "component_array[contig_num_for_kmer] %u\n", component_array[contig_num_for_kmer] );
                        //shc_log_info(shc_logname, "component_array_aux[contig_num_for_kmer] %u\n", component_array_aux[contig_num_for_kmer] );
                        comp_count[component_array[contig_num_for_kmer]]++;            
                        
                        if(metis_setup.is_multiple_partition && 
                           component_array_aux[contig_num_for_kmer] < curr_component_num
                           )                            
                            comp_count_aux[component_array_aux[contig_num_for_kmer]]++;
                        
                    }                    
                }     
                // check if this sequence
                if(!comp_count.empty())                
                {
                    auto max_iter = std::max_element(comp_count.begin(), comp_count.end(), 
                            [](const std::pair<comp_num_t, int> & p1, 
                               const std::pair<comp_num_t, int> & p2 )
                            {
                                return p1.second < p2.second;
                            }
                    );

                    comp_num_t best_comp =  max_iter->first; 
                    
                    //shc_log_info(shc_logname, "Before add to file %s \n",sequence.c_str());
                    *(files.at(best_comp)) << header << std::endl;
                    *(files.at(best_comp)) << sequence << std::endl;
                    //shc_log_info(shc_logname, "after add to file\n");                                
                    comp_count.clear();                                                    
                }                             
#ifdef LOG_CONTIG_GRAPH                
                else
                {
                    shc_log_info(shc_logname, "%s does not belong to any comp\n", header.c_str());
                }
#endif
                //for auxilary array
                if(metis_setup.is_multiple_partition && !comp_count_aux.empty())                
                {
                    auto max_iter = std::max_element(comp_count_aux.begin(), comp_count_aux.end(), 
                            [](const std::pair<comp_num_t, int> & p1, 
                               const std::pair<comp_num_t, int> & p2 )
                            {
                                return p1.second < p2.second;
                            }
                    );

                    comp_num_t best_comp =  max_iter->first; 
                    
                    //shc_log_info(shc_logname, "Before add to file %s \n",sequence.c_str());
                    *(files.at(best_comp)) << header << std::endl;
                    *(files.at(best_comp)) << sequence << std::endl;
                    //shc_log_info(shc_logname, "after add to file\n");                                
                    comp_count_aux.clear();                                                    
                }                             
#ifdef LOG_CONTIG_GRAPH                
                else
                {
                    shc_log_info(shc_logname, "%s does not belong to any comp\n", header.c_str());
                }
#endif
                sequence.clear();
            }
            header = line;   
        }
        else //line is sequence
        {     
            sequence += line;                        
        }                          
    }
    
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        (*files.at(i)).close();                
    }
    
    std::vector<std::shared_ptr<std::ofstream> >().swap(files);
    shc_log_info(shc_logname, "Finish assigning read to component\n");
}

void Contig_graph_handler::assign_kmer_to_components()
{
    shc_log_info(shc_logname, "Start assigning read to component\n");
    //creating dir
    typedef boost::filesystem::path path_t;
    typedef boost::filesystem::path dir_t;
    typedef boost::filesystem::path filename_t;
    path_t base_path = boost::filesystem::current_path();    
    std::string base_path_str(base_path.c_str()); 
    std::string dir_path = lf->output_path_str + "/components_kmer";     
    char bases[33];
    bases[kh->kmer_length] = '\0';
    
    path_t comp_kmer_path(dir_path);
    if(boost::filesystem::exists(comp_kmer_path))
        replace_directory(comp_kmer_path);
    else
        add_directory(comp_kmer_path);
            
    //create files    
    std::vector<std::shared_ptr<std::ofstream> > files;
    
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        std::string file_path = dir_path + 
                        (std::string("/kmer") + std::to_string(i));
        std::shared_ptr<std::ofstream> file(new std::ofstream);
        file->open(file_path.c_str());
        
        files.push_back(file);
    }
    shc_log_info(shc_logname, "Created kmer files\n");
    
    
    Kmer_counter_map_iterator it;
    for(it=kh->kmer_counter.begin(); it!=kh->kmer_counter.end(); it++)
    {        
        decode_kmer(bases, &(it->first), kh->kmer_length);
        *files.at(component_array[it->second.contig]) << bases << "\t" 
                << it->second.count << std::endl;   
        
        if(metis_setup.is_multiple_partition && 
        component_array_aux[it->second.contig] < curr_component_num)
        {
            *files.at(component_array_aux[it->second.contig]) << bases << "\t" 
                << it->second.count << std::endl;   
        }
    }
    
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        (*files.at(i)).close();                
    }
    std::vector<std::shared_ptr<std::ofstream> >().swap(files);
    shc_log_info(shc_logname, "Finish assigning kmer to component\n");    
}