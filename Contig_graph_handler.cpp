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
        std::ofstream kmer_outfile (component_kmer_filename_prefix + 
                               std::to_string(curr_component_num));
        
        contig_num_t root_contig = *explorable_contig_set.begin();
        contig_stack.push(root_contig);
        explorable_contig_set.erase(explorable_contig_set.begin());
        
        // add that node, so each graph has at least one node
        if(contig_vertex_map.find(root_contig) == contig_vertex_map.end())
        {
            vd_t vd = boost::add_vertex(graph);
            contig_vertex_map[root_contig] = vd;            
            graph[vd] = bundled_contig_index(root_contig);
            shc_log_info(shc_logname, "Root contig %u with vd %u\n", root_contig, vd);
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
                kmer_outfile << kmer_array << "\t" << (*it).second.count<< std::endl;
                //find new contig
                get_connected_contig_index(it);
            }
        }
        kmer_outfile.close();
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
    if( num_vertices <= metis_setup.partition_size)
    {
        vip_t vip = vertices(graph);
        //update contig to component info
        for(vi_t it=vip.first; it!=vip.second; it++)
        {
            component_array[graph[*it].contig_index] = curr_component_num;       
        }
        curr_component_num++;        
    }
    else        //call metis
    {
        create_metis_array_input();
        metis_setup.num_vertices = num_vertices;
        metis_setup.num_partition = std::min(100, 
                num_vertices/metis_setup.partition_size + 1);
        idx_t ufactor = static_cast<idx_t>(METIS_IMBALANCE_PRECISION 
                                                   * (metis_setup.overload-1));
        std::vector<idx_t> part(num_vertices);
        
        metis_setup.options[METIS_OPTION_UFACTOR] = ufactor;            
        
        int ret = METIS_PartGraphRecursive(
              &metis_setup.num_vertices, &metis_setup.ncon, 
              &metis_input.xadj.at(0), &metis_input.adjncy.at(0),
              NULL, NULL, &metis_input.adjwgt.at(0), &num_vertices, 
              NULL, NULL, metis_setup.options, &metis_setup.objval, &part.at(0));
        
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
            
            num_component_in_graph = std::max(num_component_in_graph, part[i]);                        
        }        
        num_component_in_graph = num_component_in_graph + 1;
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
        graph[vd] = bundled_contig_index(i);
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
    shc_log_info(shc_logname, "\n contig: compoent \n");
    int i=0;
    for(std::vector<contig_num_t>::iterator it=component_array.begin();
            it != component_array.end(); it++)
    {
        shc_log_info(shc_logname, "%d   %u\n", i++, *it);
    }
}

void Contig_graph_handler::run_metis_get_components()
{
    ;
}
/*
 typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits < Graph >::vertex_iterator vertex_iter;    
    
    void operator()(const Vertex& v) const
    {
        typedef boost::graph_traits<Graph> GraphTraits;
        typename boost::property_map<Graph, boost::vertex_index_t>::type 
          index = boost::get(boost::vertex_index, g);
        
        typename GraphTraits::out_edge_iterator out_i, out_end;
        typename GraphTraits::edge_descriptor e;
        
        for (boost::tie(out_i, out_end) = boost::out_edges(v, g); 
             out_i != out_end; ++out_i) 
        {
            e = *out_i;
            Vertex src = boost::source(e, g), targ = boost::target(e, g);                        
            outfile <<  index[src]+ METIS_NODE_OFFSET << " " 
                    << index[targ]+ METIS_NODE_OFFSET << " " 
                    << graph[e].weight << std::endl;            
        }        
        
        //get vertices
        std::cout << "vertices(g) = ";
        std::pair<vertex_iter, vertex_iter> vp;        
        for (vp = boost::vertices(g); vp.first != vp.second; ++vp.first)
            std::cout << index[*vp.first] << " ";
        std::cout << std::endl;
    }   
    
    Graph& g;
    std::ofstream& outfile;*/