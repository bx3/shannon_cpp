#include "Contig_graph_handler.h"

void Contig_graph_handler::group_components()
{          
    shc_log_info(shc_logname, "Start constructing contig graph\n");  
    uint8_t * contig_start = NULL;
    Kmer_counter_map_iterator it;
    char kmer_array[KMER_HOLDER_LEN];
    uint64_t byte = 0;
    kmer_array[kh->kmer_length] = '\0';
    std::string component_kmer_filename_prefix = "output/components/kmer_";    
    
    //put every contig into the set
    for(contig_num_t i=0; i<ch->num_contig; i++)
        explorable_contig_set.insert(i);
    
    while(!explorable_contig_set.empty())
    {
        shc_log_info(shc_logname, "creating kmer out file for %u\n",curr_component_num);
        //std::ofstream kmer_outfile (component_kmer_filename_prefix + 
        //                       std::to_string(curr_component_num));
        
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
            contig_start = &(ch->contig_list.at(ch->delimitor.at(i)));
            //where j is kmer index
            for(Contig_handler::size_type j=0;
                j<ch->contig_len_list.at(i)-kh->kmer_length+1;   j++)
            {
                //get each kmer and the corresponding info
                get_xmer_at_index(contig_start, ch->contig_len_list.at(i), j, 
                                                kh->kmer_length, kmer_array);                 
                encode_kmer(kmer_array, &byte, kh->kmer_length);
                it = kh->kmer_counter.find(byte);                                
                if(it == kh->kmer_counter.end())
                    shc_log_error("contain impossible kmer\n");
                
                //dump to out
                //kmer_outfile << kmer_array << "\t" << (*it).second.count<< std::endl;
                //find new contig
                get_connected_contig_index(it);
            }
        }
        //kmer_outfile.close();
        shc_log_info(shc_logname, "close kmer out file for %u\n",curr_component_num);
        //up to here we have a component and its graph
        break_and_keep_component();
        graph.clear();
        contig_vertex_map.clear();
    }
    
    shc_log_info(shc_logname, "Finish constructing contig graph\n");
}

void Contig_graph_handler::break_and_keep_component()
{
    typedef std::pair<vi_t, vi_t> vip_t;
    idx_t num_vertices = boost::num_vertices(graph);
    //shc_log_info(shc_logname, "num vertices %u, partition_size %u\n",num_vertices , metis_setup.partition_size);
    if( num_vertices <= metis_setup.partition_size)
    {
        vip_t vip = vertices(graph);
        //update contig to component info
        for(vi_t it=vip.first; it!=vip.second; it++)
        {
            component_array[graph[*it].contig_index] = curr_component_num;  
            //shc_log_info(shc_logname, "component %u has contig %u\n", curr_component_num, graph[*it].contig_index);
        }
        curr_component_num++;                
    }
    else        //call metis
    {
        shc_log_info(shc_logname, "Starting metis \n");
        create_metis_array_input();
        metis_setup.num_vertices = num_vertices;
        metis_setup.num_partition = std::min(100, 
                num_vertices/metis_setup.partition_size + 1);
        shc_log_info(shc_logname, "metis to asked to break %d partition\n", metis_setup.num_partition);
        idx_t ufactor = static_cast<idx_t>(METIS_IMBALANCE_PRECISION 
                                                   * (metis_setup.overload-1));
        std::vector<idx_t> part(num_vertices);
        
        //metis_setup.options[METIS_OPTION_UFACTOR] = ufactor;            
        
        int ret = METIS_PartGraphRecursive(
              &metis_setup.num_vertices, &metis_setup.ncon, 
              &metis_input.xadj.at(0), &metis_input.adjncy.at(0),
              NULL, NULL, &metis_input.adjwgt.at(0), &metis_setup.num_partition, 
              NULL, NULL, NULL, &metis_setup.objval, &part.at(0));
        
        if(ret == METIS_ERROR_INPUT)
            shc_log_error("METIS_ERROR_INPUT\n");
        if(ret == METIS_ERROR_MEMORY)
            shc_log_error("METIS_ERROR_MEMORY\n");
        if(ret == METIS_ERROR)
            shc_log_error("METIS_OTHER_ERROR\n");
        
        idx_t num_component_in_graph = 0;
        for(idx_t i=0; i<num_vertices; i++)
        {
            // for vertex ith descriptor, part[i] contains local component num            
            component_array[graph[static_cast<vd_t>(i)].contig_index] = 
                                                 curr_component_num + part[i];            
            shc_log_info(shc_logname, "metis component %d has contig %u\n", part[i], graph[i].contig_index);
            num_component_in_graph = std::max(num_component_in_graph, part[i]);                        
        }                
        num_component_in_graph = num_component_in_graph + 1;
        shc_log_info(shc_logname, "metis breaks to %d partition\n", num_component_in_graph);
        
        curr_component_num += num_component_in_graph;
        std::vector<idx_t>().swap(part);
    }
}

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
                    increment_edge_weight(contigA_index, contigB_index);                    
                    if (explorable_contig_set.find(contigB_index) != 
                        explorable_contig_set.end())
                    {
                        //if contigB_index is found
                        contig_stack.push(contigB_index);
                        explorable_contig_set.erase(contigB_index);
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
    if(contig_vertex_map.find(i) == contig_vertex_map.end())
    {
        vd = boost::add_vertex(graph);
        contig_vertex_map[i] = vd;   
        graph[vd] = bundled_contig_index(i);
        shc_log_info(shc_logname, "Added contig %u with vd %u\n", i, vd);
    }
    
    if(contig_vertex_map.find(j) == contig_vertex_map.end())
    {
        vd = boost::add_vertex(graph);
        contig_vertex_map[j] = vd;        
        graph[vd] = bundled_contig_index(j);
        shc_log_info(shc_logname, "Added contig %u with vd %u\n", j, vd);
    }
        
    std::pair<ed_t, bool> aer;
        
    aer = boost::edge(contig_vertex_map[i], contig_vertex_map[j], graph);    
    if(aer.second)
    {                
        graph[aer.first].weight++;
        shc_log_info(shc_logname, "edge between %u %u : weight increased to %u\n",
                i, j ,graph[aer.first].weight);        
    }
    else //add edge
    {
        aer = boost::add_edge(contig_vertex_map[i], contig_vertex_map[j], graph);
        assert(aer.second);
        shc_log_info(shc_logname, "Edge between %u %u\n", i, j);
        graph[aer.first] = bundled_weight(1);                
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
}

void Contig_graph_handler::create_metis_format_from_graph()
{
    shc_log_info(shc_logname, "Start making METIS format file\n");    
    typename boost::graph_traits < graph_t >::adjacency_iterator vi, vi_end;    
    vd_t vd;
    std::ofstream outfile (metis_filename);
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
    
    log_metis_input_data();
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
    info_log_str_without_new_line(shc_logname,"\n");
}

void Contig_graph_handler::log_component_array()
{
    shc_log_info(shc_logname, "contig: %u component \n", curr_component_num);
    int i=0;
    for(std::vector<contig_num_t>::iterator it=component_array.begin();
            it != component_array.end(); it++)
    {
        shc_log_info(shc_logname, "%d   %u\n", i++, *it);
    }
    shc_log_info(shc_logname, "Finish log_component_array\n");
}


void Contig_graph_handler::assign_reads_to_components(
                                std::string& read_filename, int num_test )
{
    //Fasta_reader fasta_reader("SE_read.fasta");        
    shc_log_info(shc_logname, "Start assigning read to component\n");
    //creating dir
    typedef boost::filesystem::path path_t;
    typedef boost::filesystem::path dir_t;
    typedef boost::filesystem::path filename_t;
    path_t base_path = boost::filesystem::current_path();    
    dir_t comp_read_local("output/components_reads");
    path_t dir_path = base_path / comp_read_local;  
    
    add_directory(dir_path);
    shc_log_info(shc_logname, "Created dir %s\n", dir_path.c_str());
            
    //create files    
    std::vector<std::shared_ptr<std::ofstream> > files;
    
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        path_t file_path = dir_path / 
                        (std::string("comp") + std::to_string(i));
        std::shared_ptr<std::ofstream> file(new std::ofstream);
        file->open(file_path.c_str());
        
        files.push_back(file);
    }
    shc_log_info(shc_logname, "Created contig files\n");
    
    //read and process files
    std::ifstream file_reader(read_filename);  
    shc_log_info(shc_logname, "successfully open read file\n");
    std::string line;
    std::string header;
    std::string sequence;    
    int interval = 1;
    uint64_t byte = 0;
    
    std::map<comp_num_t, int> comp_count;
        
    while (file_reader.good() ) 
    {
        std::getline(file_reader, line);
        if(line[0] == '>')
        {
            header = line;
            if(!sequence.empty())
            {
                //process
                
                interval = (sequence.size()-kh->kmer_length+1)/num_test;
                for(int i=0; i<num_test; i++)
                {
                    //shc_log_info(shc_logname, "Before encoding readsize %u at %d\n",sequence.size(), i*interval);
                    encode_kmer(&sequence.at(0)+i*interval, &byte, kh->kmer_length);
                    //shc_log_info(shc_logname, "After encoding\n");
                    comp_num_t comp_num = 
                                component_array[kh->kmer_counter[byte].contig];
                    comp_count[comp_num]++;            
                }        
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
                sequence.clear();
            }
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
    dir_t comp_read_local("output/components_kmer");
    path_t dir_path = base_path / comp_read_local;  
    char bases[33];
    bases[kh->kmer_length] = '\0';
    
    add_directory(dir_path);
    shc_log_info(shc_logname, "Created dir %s\n", dir_path.c_str());
            
    //create files    
    std::vector<std::shared_ptr<std::ofstream> > files;
    
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        path_t file_path = dir_path / 
                        (std::string("kmer") + std::to_string(i));
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
    }
    
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        (*files.at(i)).close();                
    }
    std::vector<std::shared_ptr<std::ofstream> >().swap(files);
    shc_log_info(shc_logname, "Finish assigning kmer to component\n");    
}
