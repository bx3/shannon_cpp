#include "Sparse_flow_handler.h"

double norm_1_helper(double x, double y) {return x + abs(y);}

template<class T1, class T2> double dot_product(T1 & v1, T2 & v2)
{
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

void * Sparse_flow_handler::run_sparse_flow_handler_helper( Comp_graph comp_graph,
            std::string node_path, std::string edge_path, std::string path_path,
            pthread_mutex_t * write_lock_ptr, std::string output_path)
{
    comp_id = comp_graph.comp_i;
    graph_id = comp_graph.graph_i;
    if(!load_graph(node_path, edge_path))
    {
        clear();
        return NULL;
    }

    load_path(path_path);

    //std::cout << "get_num_nodes " << get_num_nodes() << std::endl;
    //std::cout << "get_num_nodes " << get_num_edges() << std::endl;

    add_start_end_node();

    if (get_num_nodes() <= 3)
    {
        dump_SF_seq(output_path, write_lock_ptr);
        clear();
    }
    else
    {
        //log_graph_struct(false);
        sparse_flow_on_all_nodes();
        dump_SF_seq(output_path, write_lock_ptr);
        clear();
    }
}


void Sparse_flow_handler::clear()
{
    graph.clear();
    sorted_vds.clear();
    sorted_vds.shrink_to_fit();
    nodeId_vd_map.clear();
    known_paths.clear();
    known_paths.shrink_to_fit();
    paths_for_node.clear();
    vd_contains_map.clear();

    vd_start = NULL;
    vd_end = NULL;
    glob_node_id = 0;
}

// process one component
void Sparse_flow_handler::process_one_graph_component(int graph_i,std::string & node_path,
                        std::string & edge_path, std::string & path_path)
{
    std::cout << "process_one_graph_component, not multithread" << std::endl;
    if(!load_graph(node_path, edge_path))
    {
        shc_log_info(shc_logname, "invalid input graph\n", graph_i);
        return;
    }
#ifdef LOG_SF
    std::cout << "after load graph, num nodes " << get_num_nodes() << " num edges "
              << get_num_edges() << std::endl;
#endif

    comp_id = 3000;
    graph_id = graph_i;

    //exit(0);

    //log_graph_struct(false);
    load_path(path_path);

    //std::cout << "get_num_nodes " << get_num_nodes() << std::endl;
    //std::cout << "get_num_nodes " << get_num_edges() << std::endl;

    add_start_end_node();
#ifdef LOG_SF
    std::cout << "after add start end, node num"<< get_num_nodes() << " num edges "
              << get_num_edges() << std::endl;
#endif
    //log_graph_struct(false);
    if(get_num_nodes() <= 3)
    {
        dump_SF_seq(setting.local_files.reconstructed_seq_path);
        return;
    }

    //log_all_node(true);

    sparse_flow_on_all_nodes();
#ifdef LOG_SF
    std::cout << "after SF, node num " << get_num_nodes() << " num edges "
              << get_num_edges() << std::endl;
    log_all_node(true);
#endif
    //log_graph_struct(false);
    /*
    if(!acyclic_check())
    {
        shc_log_error("after sparse flow, graph becomes cyclic\n");
        exit(1);
    }
    else
    {
        shc_log_info(shc_logname, "after sparse flow, graph is acyclic\n");
    }
    */
    num_out_seq = 0;
    dump_SF_seq(setting.local_files.reconstructed_seq_path);
#ifdef LOG_SF
    shc_log_info(shc_logname, "Finish processing graph components %d\n", graph_i);
#endif
}

// run everything in one thread
void Sparse_flow_handler::process_all_graph(int comp_i)
{
#ifdef LOG_SF
    shc_log_info(shc_logname, "Start processing all seq graph components\n");
#endif
    int i = 0;
    Local_files & lf = setting.local_files;

    std::string node_path_prefix(lf.output_seq_graph_path +
                                 lf.node_prefix + std::to_string(comp_i) +
                                 "/node");
    std::string edge_path_prefix(lf.output_seq_graph_path +
                                 lf.edge_prefix + std::to_string(comp_i) +
                                 "/edge");
    std::string path_path_prefix(lf.output_seq_graph_path +
                                 lf.path_prefix + std::to_string(comp_i) +
                                 "/path");
    bool flag = boost::filesystem::exists(node_path_prefix + std::to_string(i)) &&
                boost::filesystem::exists(edge_path_prefix + std::to_string(i)) &&
                boost::filesystem::exists(path_path_prefix + std::to_string(i));
    while(flag)
    {
        std::string graph_node_path = lf.node_prefix + std::to_string(i);
        std::string graph_edge_path = lf.edge_prefix + std::to_string(i);
        std::string graph_path_path = lf.path_prefix + std::to_string(i);

        process_one_graph_component(i, graph_node_path, graph_edge_path, graph_path_path);

        i++;
        flag = boost::filesystem::exists(node_path_prefix + std::to_string(i)) &&
               boost::filesystem::exists(edge_path_prefix + std::to_string(i)) &&
               boost::filesystem::exists(path_path_prefix + std::to_string(i));
    }
#ifdef LOG_SF
    shc_log_info(shc_logname, "Finish processing %d seq graph components\n", i);
#endif
}

//bool load_py_graph(std::string & node_file, std::string & edge_file)


bool Sparse_flow_handler::load_graph(std::string & node_file, std::string & edge_file)
{
#ifdef LOG_SF
    shc_log_info(shc_logname, "Start load graph files\n");
#endif
    if(!boost::filesystem::exists(node_file))
    {
        std::cout << "input comp " << comp_id << " graph " << graph_id
                  <<  "node file does not exist" << std::endl;
        std::cout << node_file << std::endl;
        return false;
    }
    if(!boost::filesystem::exists(edge_file))
    {
        std::cout << "input graph edge file does not exist" << std::endl;
        std::cout << edge_file << std::endl;
        return false;
    }

    std::ifstream node_file_reader(node_file);
    std::string header_str;
    std::getline(node_file_reader, header_str);
    std::vector<std::string> strs;
    boost::split(strs, header_str, boost::is_any_of("\t"));
    node_file_reader.close();

    bool is_valid_node = false;
    bool is_valid_edge = false;

    if(strs.size()== 3)
    {
        is_valid_node = load_node(node_file);
        is_valid_edge = load_edge(edge_file);
    }
    else if(strs.size()== 4)
    {
        is_valid_node = load_py_node(node_file);
        is_valid_edge = load_py_edge(edge_file);
    }
    else
    {
        shc_log_error("comp %d graph %d node file header read error: %s\n",
                        comp_id, graph_id, header_str.c_str());
        return false;
    }

    if(!is_valid_node)
    {
        std::cout << "input graph is not a valid graph, nodes file error" << std::endl;
    }
    if(!is_valid_edge)
    {
        std::cout << "input graph is not a valid graph, edge file error" << std::endl;
    }

    //std::cout << "Graph contains " << (boost::num_vertices(graph)) << " nodes" << std::endl;
    //std::cout << "Graph contains " << (boost::num_edges(graph)) << " edges" << std::endl;

#ifdef LOG_SF
    shc_log_info(shc_logname, "node file %s\n", node_file.c_str());
    shc_log_info(shc_logname, "edge file %s\n", edge_file.c_str());

    shc_log_info(shc_logname, "Graph: %d nodes, %d edges\n",
                            boost::num_vertices(graph), boost::num_edges(graph));
    shc_log_info(shc_logname, "Finish load graph files\n");
#endif
    return is_valid_node && is_valid_edge;
}

/**
*  Nodes are assumed to be loaded in topologocal order
*/
bool Sparse_flow_handler::load_node(std::string & node_file)
{
    std::ifstream file_reader(node_file);
    std::string node_id_str, bases, count_str;
    node_id_t node_id;
    kmer_count_node_t count;

    int num_node = 0;

    std::getline(file_reader, node_id_str); // get first line
    while(  std::getline(file_reader, node_id_str, '\t') &&
            std::getline(file_reader, bases, '\t') &&
            std::getline(file_reader, count_str))
    {
        node_id = boost::lexical_cast<node_id_t>(node_id_str);
        count = boost::lexical_cast<kmer_count_node_t>(count_str);
        vd_t vd = boost::add_vertex(graph);
        sorted_vds.push_back(vd);
        //shc_log_info(shc_logname, "VD: %d\n", node_id);
        graph[vd] = bundled_node_p(bases, node_id);
        nodeId_vd_map.insert(std::make_pair(node_id, vd));
        if (glob_node_id < node_id)
            glob_node_id = node_id;

        num_node++;
    }
    glob_node_id++;
    //std::cout << "glob_node_id is " <<glob_node_id << std::endl;
    file_reader.close();
    return num_node > 0;
}

bool Sparse_flow_handler::load_py_node(std::string & node_file)
{
    std::ifstream file_reader(node_file);
    std::string node_id_str, bases, count_str, norm_str;
    node_id_t node_id;
    kmer_count_node_t count;
    uint64_t norm;

    int num_node = 0;

    std::getline(file_reader, node_id_str); // get first line
    while(  std::getline(file_reader, node_id_str, '\t') &&
            std::getline(file_reader, bases, '\t') &&
            std::getline(file_reader, count_str, '\t') &&
            std::getline(file_reader, norm_str) )
    {
        node_id = boost::lexical_cast<node_id_t>(node_id_str);
        count = boost::lexical_cast<double>(count_str);
        vd_t vd = boost::add_vertex(graph);
        sorted_vds.push_back(vd);
        //shc_log_info(shc_logname, "VD: %d\n", node_id);
        graph[vd] = bundled_node_p(bases, node_id);
        nodeId_vd_map.insert(std::make_pair(node_id, vd));
        if (glob_node_id < node_id)
            glob_node_id = node_id;

        num_node++;
    }
    glob_node_id++;
    //std::cout << "glob_node_id is " <<glob_node_id << std::endl;
    file_reader.close();
    return num_node > 0;
}

bool Sparse_flow_handler::load_edge(std::string & edge_file)
{
    std::ifstream file_reader(edge_file);
    std::string node_from_str, node_to_str, weight_str, count_str;
    node_id_t from_node_id, to_node_id;
    double count;
    edge_weight_t weight;

    std::getline(file_reader, node_from_str); // get first line

    int num_edge = 0;

    while(  std::getline(file_reader, node_from_str, '\t') &&
            std::getline(file_reader, node_to_str, '\t') &&
            std::getline(file_reader, weight_str, '\t') &&
            std::getline(file_reader, count_str))
    {
        from_node_id = boost::lexical_cast<node_id_t>(node_from_str);
        to_node_id = boost::lexical_cast<node_id_t>(node_to_str);
        weight = boost::lexical_cast<edge_weight_t>(weight_str);
        count = boost::lexical_cast<double>(count_str);
        vd_t vd_from = nodeId_vd_map[from_node_id];
        vd_t vd_to = nodeId_vd_map[to_node_id];
        aer_t aer = boost::add_edge(vd_from, vd_to, graph);
        graph[aer.first] = bundled_edge_p(weight, count);
        num_edge++;
        //shc_log_info(shc_logname, "From VD: %d TO VD: %d\n", from_node_id, to_node_id);
    }
    file_reader.close();
    return num_edge > 0 ;
}

bool Sparse_flow_handler::load_py_edge(std::string & edge_file)
{
    std::ifstream file_reader(edge_file);
    std::string node_from_str, node_to_str, weight_str, count_str, norm_str;
    node_id_t from_node_id, to_node_id;
    double count;
    edge_weight_t weight;

    std::getline(file_reader, node_from_str); // get first line

    int num_edge = 0;

    while(  std::getline(file_reader, node_from_str, '\t') &&
            std::getline(file_reader, node_to_str, '\t') &&
            std::getline(file_reader, weight_str, '\t') &&
            std::getline(file_reader, count_str, '\t') &&
            std::getline(file_reader, norm_str) )
    {
        from_node_id = boost::lexical_cast<node_id_t>(node_from_str);
        to_node_id = boost::lexical_cast<node_id_t>(node_to_str);
        weight = boost::lexical_cast<edge_weight_t>(weight_str);
        count = boost::lexical_cast<double>(count_str);
        vd_t vd_from = nodeId_vd_map[from_node_id];
        vd_t vd_to = nodeId_vd_map[to_node_id];
        aer_t aer = boost::add_edge(vd_from, vd_to, graph);
        graph[aer.first] = bundled_edge_p(weight, count);
        num_edge++;
        //shc_log_info(shc_logname, "From VD: %d TO VD: %d\n", from_node_id, to_node_id);
    }
    file_reader.close();
    return num_edge > 0 ;
}

void Sparse_flow_handler::load_path(std::string & path_file)
{
#ifdef LOG_SF
    shc_log_info(shc_logname, "load_path\n");
#endif
    std::ifstream file_reader(path_file);
    std::string line;
    std::getline(file_reader, line);
    node_id_t node_id;
    int i = 0;
    while(std::getline(file_reader, line) )
    {
        //shc_log_info(shc_logname, "new line\n");
        std::vector<std::string> tokens;
        std::vector<vd_t> nodes_in_path;
        boost::algorithm::split(tokens, line, boost::is_any_of("\t"));
        //for(int i=0; i< tokens.size(); i++)
        //{
        //    shc_log_info(shc_logname, "print string %d: %s\n", i, tokens[i].c_str());
        //}

        int location = 0;
        for(std::string & node_id_str: tokens)
        {
            if(node_id_str.empty() || node_id_str[0]==' ')
                continue;
            //shc_log_info(shc_logname, "%s\n", node_id_str.c_str());
            node_id = boost::lexical_cast<node_id_t>(node_id_str);
            //shc_log_info(shc_logname, "node id %d\n", node_id);
            vd_t vd = nodeId_vd_map[node_id];
            nodes_in_path.push_back(vd);
            Node_path_MMap_iterator_pair np_it_pair =
                    paths_for_node.equal_range(vd);
            if(np_it_pair.first == np_it_pair.second) // no such entry
            {
                //shc_log_info(shc_logname, "path for node location %d\n", location);
                paths_for_node.insert(std::make_pair(vd, Path_info(i, location)));
            }
            else
            {
                int diff = std::distance(np_it_pair.first, np_it_pair.second);
                //shc_log_info(shc_logname, "diff %d\n", diff);
                if(diff < sf_setting.path_sparsity)
                {
                //    shc_log_info(shc_logname, "path for node location %d\n", location);
                    paths_for_node.insert(std::make_pair(vd, Path_info(i, location)));
                }
            }
            location++;
        }
        i++;
        known_paths.push_back(nodes_in_path);
    }
    file_reader.close();
}

void Sparse_flow_handler::add_start_end_node()
{
#ifdef LOG_SF
    shc_log_info(shc_logname, "Start adding start,end nodes\n");
#endif
    std::vector<vd_t> front_vd, back_vd;
    vip_t vip = boost::vertices(graph);
    for(vi_t vi=vip.first; vi!=vip.second; vi++)
    {
        vd_t vd = *vi;
        if(boost::in_degree(vd, graph) == 0)
            front_vd.push_back(vd);
        if(boost::out_degree(vd, graph) == 0)
            back_vd.push_back(vd);
    }

    vd_start = boost::add_vertex(graph);
    vd_end = boost::add_vertex(graph);
    graph[vd_start] = bundled_node_p("Start_", NODE_ID_START);// NODE_ID_START
    graph[vd_end] = bundled_node_p("End_", NODE_ID_END);// NODE_ID_END

    for(std::vector<vd_t>::iterator it =front_vd.begin();
                                    it!=front_vd.end(); it++)
    {
        bundled_node_p & curr_node = graph[*it];
        aer_t aer = boost::add_edge(vd_start, *it, graph);
        graph[aer.first] = bundled_edge_p(0, 0);
    }

    for(std::vector<vd_t>::iterator it =back_vd.begin();
                                    it!=back_vd.end(); it++)
    {
        bundled_node_p & curr_node = graph[*it];
        aer_t aer = boost::add_edge(*it, vd_end, graph);
        graph[aer.first] = bundled_edge_p(0, 0);
    }
#ifdef LOG_SF
    shc_log_info(shc_logname, "Finish adding start,end nodes\n");
#endif
}

void Sparse_flow_handler::sparse_flow_on_all_nodes()
{
#ifdef LOG_SF
    shc_log_info(shc_logname, "\t\tStart SF\n");
#endif
    bool done = false;
    bool is_concern_unqiue = false;
    // all nodes contains itself
    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        std::deque<vd_t> contain_array(1, vd);
        vd_contains_map.insert(std::make_pair(vd, contain_array));
    }

    std::ofstream writer("/home/bowen/Shannon_D_RNA_seq/sf.log");

    int algo_iter = 0;
    int node_iter = 0;
    while(!done)
    {
#ifdef LOG_SF
        shc_log_info(shc_logname, "new algo iteration\n");
#endif
        //std::cout << "algo iter" << algo_iter
        //          << " num nodes " << sorted_vds.size() << std::endl;
        // modify later in topological order
        std::set<vd_t> remove_vd;

        //link_all_Y_path(remove_vd);

        //if (algo_iter ==1)
        //{
        //    std::cout << "node_iter " << node_iter << std::endl;
        //    std::cout << "node left " << sorted_vds.size() << std::endl;
        //    exit(0);
        //}
        node_iter =  0;

        for(uint64_t node_i=0; node_i<sorted_vds.size(); node_i++)
        {
            vd_t vd = sorted_vds[node_i];
            //shc_log_info(shc_logname, "***********************\n");
            //log_local_node(vd, false);
            //shc_log_info(shc_logname, "anS %d %s\n", graph[vd].node_id, graph[vd].seq.c_str());
            //continue;
            //node_iter++;
            if(vd==vd_start || vd==vd_end)
            {
                //shc_log_info(shc_logname, "node is source or sink\n");
                continue;
            }
            if(boost::in_degree(vd, graph) <=1)
            {
                //shc_log_info(shc_logname, "in degree <= 1\n");
                continue;
            }
            //if(boost::out_degree(vd, graph) ==1)
            //{
            //    out_eip_t out_eip = boost::out_edges(vd,graph);
            //    vd_t vd_target = boost::target(*(out_eip.first), graph);
            //    if (vd_target == vd_end)
            //        continue;
            //}

            //node_i++;

            if(boost::out_degree(vd, graph) ==1)
            {
                //shc_log_info(shc_logname, "Y node modify\n");
                std::vector<vd_t> surround_nodes;
                out_eip_t out_eip = boost::out_edges(vd,graph);
                vd_t vd_target = boost::target(*out_eip.first, graph);
                in_eip_t in_eip = boost::in_edges(vd,graph);
                for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
                {
                    vd_t vd_source = boost::source(*in_ei, graph);
                    //if (vd_source != vd_start)
                    //{
                    //    shc_log_info(shc_logname, "%d -> %d\n",
                    //            graph[vd_source].node_id, graph[vd_target].node_id);
                    //}
                }


                link_Y_path(vd);
                //condense_nodes(surround_nodes, remove_vd);
            }
            else
            {

                //shc_log_info(shc_logname, "X node modify\n");
                node_iter++;
                std::vector<double> flow;
                std::vector<vd_t> in_nodes, out_nodes;
                std::vector<double> in_counts, out_counts;
                get_in_out_edges(vd, in_nodes, out_nodes, in_counts, out_counts);

                // get known path indicator in an array
                 // for indicating which sub flow is unknown
                std::vector<int> sub_flow_indicator;
                get_sub_flow_indicator(vd, sub_flow_indicator, in_nodes, out_nodes);


                //for(double in_c : in_counts)
                //    printf("in count :  %f\n", in_c);
                //for(double out_c : out_counts)
                //    printf("out: %f\n", out_c);

                //exit(0);
                //shc_log_info(shc_logname, "BSF\n");
                //log_local_node(vd, false);
                bool is_unique = sparse_flow_on_one_node(vd, flow,
                                   in_counts, out_counts, sub_flow_indicator, in_nodes, out_nodes);
                //if(!is_unique)
                //    shc_log_warning("non unique flow\n");
#ifdef LOG_SF
                int num_out = boost::out_degree(vd, graph);
                for(int k=0; k<flow.size(); ++k)
                {
                    int i, j;
                    get_i_j(k, num_out, i, j);
                    if(flow[k]!=0)
                    {
                        vd_t u = in_nodes[i];
                        vd_t w = out_nodes[j];
                        //std::cout << "flow from " << graph[u].node_id << " to "
                        //          << graph[w].node_id << std::endl;
                        shc_log_info(shc_logname, "***********************\n");
                        shc_log_info(shc_logname, "%d -> %d\n", graph[u].node_id, graph[w].node_id);
                    }
                }
                shc_log_info(shc_logname, "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n");
#endif

                //modify_writer << "vd " << graph[vd].node_id
                //          << " num_in: " << in_counts.size()
                //          << ", num_out: " << out_counts.size()
                //          << " non zero flow " << get_norm_0(flow)
                //          << std::endl;
                modify_local_graph(vd, flow, in_nodes, out_nodes);

                //link_leaf_node_to_terminal(in_nodes, out_nodes);
                //condense_nodes(in_nodes,  remove_vd);
                //condense_nodes(out_nodes, remove_vd);
                remove_vd.insert(vd);
                //exit(0); //hello
            }
        }
        //log_term_node(true, true, false);
        //file_node << "node iter" << std::to_string(node_iter) << std::endl;
        //std::cout << "SF " << node_iter << " nodes" << std::endl;
        algo_iter++;

        //condense_graph(remove_vd);

        //if (sf_setting.is_use_topo_reverse)
        //    update_sorted_vds(true);
        //else
            update_sorted_vds(false);

        remove_decomposed_node(remove_vd);

#ifdef LOG_SF
        shc_log_info(shc_logname, "\t\t%d SF run, %d node remain to process\n",
                                                algo_iter, sorted_vds.size());
#endif
        //exit(0);
        if(sorted_vds.size() == 0)
            done = true;
    }

#ifdef LOG_SF
    shc_log_info(shc_logname, "Finish SF\n");
#endif
    //check_any_condensible_node();
    //log_graph_struct(false);

    //file_node.close();
}

void Sparse_flow_handler::log_flow(std::vector<double> & flow,
                std::vector<vd_t> & in_nodes, std::vector<vd_t> & out_nodes)
{
    int num_out = out_nodes.size();
    for(int k=0; k<flow.size(); ++k)
    {
        int i, j;
        get_i_j(k, num_out, i, j);
        if(flow[k]!=0)
        {
            vd_t u = in_nodes[i];
            vd_t w = out_nodes[j];
            //std::cout << "flow from " << graph[u].node_id << " to "
            //          << graph[w].node_id << std::endl;
            shc_log_info(shc_logname, "***********************\n");
            shc_log_info(shc_logname, "%d -> %d\n", graph[u].node_id, graph[w].node_id);
        }
    }
    shc_log_info(shc_logname, "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n");
}

void Sparse_flow_handler::link_all_Y_path(std::set<vd_t> & remove_vd)
{
    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
#ifdef LOG_SF
        shc_log_info(shc_logname, "*********************************\n");
        shc_log_info(shc_logname, "look at node %d\n", graph[vd].node_id);
        log_node_info(vd, true);
#endif
        if(vd==vd_start || vd==vd_end)
        {
#ifdef LOG_SF
            shc_log_info(shc_logname, "is source or sink\n");
#endif
            continue;
        }
        if(boost::in_degree(vd, graph)==0 && boost::out_degree(vd, graph)==0)
        {
#ifdef LOG_SF
            shc_log_info(shc_logname, "isolated node\n");
#endif
            continue;
        }

        if(local_condense(vd, remove_vd)!=NULL)
        {
#ifdef LOG_SF
            shc_log_info(shc_logname, "condensed\n");
            shc_log_info(shc_logname, "look at node %d\n", graph[vd].node_id);
            log_node_info(vd, true);
#endif
        }
        else
        {
#ifdef LOG_SF
            shc_log_info(shc_logname, "not condensed\n");
#endif
        }

        int in_degree = boost::in_degree(vd, graph);
        int out_degree = boost::out_degree(vd, graph);

        if(in_degree==1 && out_degree==1)
        {
#ifdef LOG_SF
            shc_log_info(shc_logname, "in out 1\n");
#endif
        }
        else if(in_degree==1 || out_degree==1 )
        {
#ifdef LOG_SF
            vd_t vd_front=NULL;
            vd_t vd_behind=NULL;
            if(in_degree==1)
            {
                in_eip_t in_eip = boost::in_edges(vd,graph);
                vd_front = boost::source(*in_eip.first, graph);
            }
            else
            {
                out_eip_t out_eip = boost::out_edges(vd,graph);
                vd_behind = boost::target(*out_eip.first, graph);
            }

#endif

            if(link_Y_path(vd))
            {
#ifdef LOG_SF
                shc_log_info(shc_logname, "Y linked\n");
                if(vd_front!=NULL)
                    log_node_info(vd_front, true);
                else
                    log_node_info(vd_behind, true);
#endif
                remove_vd.insert(vd);
            }
            else
            {
#ifdef LOG_SF
                shc_log_info(shc_logname, "connected to source or sink\n");
#endif
            }
        }
    }
}

void Sparse_flow_handler::
condense_nodes(std::vector<vd_t> & nodes, std::set<vd_t> & remove_vd)
{
    for(std::vector<vd_t>::iterator it= nodes.begin();
                                    it!=nodes.end(); it++)
    {
        local_condense(*it, remove_vd);
    }
}

double Sparse_flow_handler::
get_scale(std::vector<double> & in_counts, std::vector<double> & out_counts)
{
    double in_max = *std::max_element(in_counts.begin(), in_counts.end());
    double out_max = *std::max_element(out_counts.begin(), out_counts.end());
    return std::max(std::max(in_max, out_max), 0.000000000001) * 0.01;
}

bool Sparse_flow_handler::
sparse_flow_on_one_node(vd_t vd, std::vector<double> & min_flow,
              std::vector<double> & in_counts, std::vector<double> & out_counts,
              std::vector<int> & sub_flow_indicator, std::vector<vd_t> & in_nodes,
                                                        std::vector<vd_t> & out_nodes)
{
#ifdef LOG_SF
    shc_log_info(shc_logname, "\t\t SF on node %d\n", graph[vd].node_id);
#endif

    int num_in = boost::in_degree(vd, graph);
    int num_out = boost::out_degree(vd, graph);

    if(num_in == 1)
    {
        //shc_log_error("look at a  node with in edge num = 1\n");
        //exit(0);
    }
    //if(num_out == 1)
    //{
    //
    //}

    //setting
    sf_setting.tol = get_norm_1(in_counts)*sf_setting.eps;
    set_trials(num_in, num_out);

#ifdef LOG_LP_SUMMARY
    shc_log_info(shc_logname, "in nodes order:\n");
    for(int i=0; i<in_nodes.size(); i++)
        shc_log_info(shc_logname, "%u:\n", graph[in_nodes[i]].node_id);
    shc_log_info(shc_logname, "out nodes order:\n");
    for(int i=0; i<out_nodes.size(); i++)
        shc_log_info(shc_logname, "%u:\n", graph[out_nodes[i]].node_id);

    shc_log_info(shc_logname, "BEFORE in counts:\n");
    for(std::vector<double>::iterator it=in_counts.begin(); it!=in_counts.end(); ++it)
        shc_log_info(shc_logname, "%f\n", *it);
    shc_log_info(shc_logname, "BEFORE out counts:\n");
    for(std::vector<double>::iterator it=out_counts.begin(); it!=out_counts.end(); ++it)
        shc_log_info(shc_logname, "%f\n", *it);
#endif

    /*
    modify_writer << "num_in: " << num_in <<  ", num_out: " << num_out << std::endl;
    modify_writer << "in counts" << std::endl;
    for(std::vector<double>::iterator it=in_counts.begin(); it!=in_counts.end(); ++it)
        modify_writer <<  *it << std::endl;

    modify_writer << "out counts:" << std::endl;
    for(std::vector<double>::iterator it=out_counts.begin(); it!=out_counts.end(); ++it)
        modify_writer << *it << std::endl;
    modify_writer << "in sum " << sum_vector(in_counts)
                  << ", out sum " << sum_vector(out_counts) << std::endl ;
    */
    // make flow equal, and construct equal constraint
    make_in_out_flow_equal(in_counts, out_counts);
    std::vector<double> equ_constraint(in_counts);
    equ_constraint.insert(equ_constraint.end(),
                          out_counts.begin(), out_counts.end());


    //modify_writer << "after net flow" << std::endl;
    //modify_writer << "num_in: " << num_in <<  ", num_out: " << num_out << std::endl;
    //modify_writer << "in counts" << std::endl;
    //for(std::vector<double>::iterator it=in_counts.begin(); it!=in_counts.end(); ++it)
    //  modify_writer <<  *it << std::endl;

    //modify_writer << "out counts:" << std::endl;
    //for(std::vector<double>::iterator it=out_counts.begin(); it!=out_counts.end(); ++it)
    //  modify_writer << *it << std::endl;
    //modify_writer << "in sum " << sum_vector(in_counts)
    //            << ", out sum " << sum_vector(out_counts) << std::endl ;

#ifdef LOG_LP_SUMMARY
    shc_log_info(shc_logname, "num_in: %d, num_out: %d\n", num_in, num_out);
    shc_log_info(shc_logname, "in counts:\n");
    for(std::vector<double>::iterator it=in_counts.begin(); it!=in_counts.end(); ++it)
        shc_log_info(shc_logname, "%f\n", *it);
    shc_log_info(shc_logname, "out counts:\n");
    for(std::vector<double>::iterator it=out_counts.begin(); it!=out_counts.end(); ++it)
        shc_log_info(shc_logname, "%f\n", *it);
    shc_log_info(shc_logname, "in sum %f, out sum %f\n", sum_vector(in_counts), sum_vector(out_counts));
#endif



    int min_num_flow = num_in*num_out + 1;
    int flow_times = 0;
    double min_unknown_flow = 0;

#ifdef LOG_SF
    shc_log_info(shc_logname, "Before path decom %d trails, %d in %d out\n",
                                        sf_setting.num_trials, num_in, num_out);
#endif

    double scale = get_scale(in_counts, out_counts);
    //std::cout << "scale " << scale << std::endl;
    for(std::vector<double>::iterator it=equ_constraint.begin();
                            it!=equ_constraint.end(); it++)
    {
        (*it) /= scale;
    }

    // loop over trial times and choose the best
    for(int i=0; i< sf_setting.num_trials; ++i)
    {

        std::vector<double> gamma;
        std::vector<double> flow;
        // get random gamma
        for(int j=0; j<num_in*num_out; j++)
        {
            gamma.push_back(std::abs(distribution(generator)));
            //shc_log_info(shc_logname, "g %f\n", gamma.back());
        }

        double obj_value = path_decompose(vd, gamma, equ_constraint,
                                           sub_flow_indicator, flow, scale);




        //shc_log_info(shc_logname, "flow size %d\n", flow.size());
        std::vector<double> flow_sol(flow); // a copy

        //truncate_flow_value(num_in, num_out, flow, flow_sol,
        //                                        in_counts, out_counts);

        old_truncate_flow_value(num_in, num_out, flow, flow_sol,
                                                in_counts, out_counts);
#ifdef LOG_LP_SUMMARY
        shc_log_info(shc_logname, "after trunc solution summary\n");
        for(int k=0; k<flow.size(); k++)
        {
            int i, j;
            get_i_j(k, num_out, i, j);
            shc_log_info(shc_logname, "f(%d,%d) = %f\n", i,j, flow.at(k));
        }
#endif


        //modify_writer << "after trunc solution summary " << i << std::endl;
        //for(int k=0; k<flow.size(); k++)
        //{
        //    int r, j;
        //    get_i_j(k, num_out, r, j);

            //<< std::setprecision(2)
        //    modify_writer << "f(" << r << "," << j << ") = " << flow.at(k) << std::endl;
        //}

        // update procedures
        int num_flow = get_num_nonzero_flow(flow_sol, sub_flow_indicator);
        //shc_log_info(shc_logname, "num non known_path flow %d\n", num_flow);
        //modify_writer << "num_flow " << num_flow << " min_num_flow " << min_num_flow << std::endl;
        if(num_flow < min_num_flow)
        {
            //modify_writer << "\t\tupdate less "  << graph[vd].node_id<< std::endl;
            min_num_flow = num_flow;
            min_flow = flow;
            //log_flow(flow, in_nodes, out_nodes);
            flow_times = 1;
            min_unknown_flow = dot_product(sub_flow_indicator, min_flow);
        }
        else if(num_flow == min_num_flow)
        {
            std::vector<double> diff_result;
            std::transform(min_flow.begin(), min_flow.end(), flow.begin(),
                               std::back_inserter(diff_result),
                               std::minus<double>());
            if(get_norm_2(diff_result) > sf_setting.tol)
                flow_times++;

            //modify_writer << "min_unknown_flow " << min_unknown_flow << std::endl;
            if(is_update_flow(flow, min_flow, sub_flow_indicator, min_unknown_flow))
            {
                //modify_writer << "\t\tupdate equal " << graph[vd].node_id << std::endl;
                min_flow = flow;
                //log_flow(flow, in_nodes, out_nodes);
                min_unknown_flow = dot_product(sub_flow_indicator, min_flow);
            }
        }
    }


    //modify_writer << "middle check  vd: " << graph[vd].node_id
    //              << " in " << num_in << " out " << num_out << std::endl;
    //for(int k=0; k<min_flow.size(); k++)
    //{
    //    int i, j;
    //    get_i_j(k, num_out, i, j);

        //<< std::setprecision(2)
    //    modify_writer << "f(" << i << "," << j << ") = " << min_flow.at(k) << std::endl;
    //}


    //int num_flow_before_sparser = get_num_nonzero_flow(min_flow);
    //std::cout << "before py_make_flow_sparser, num non-zero flow " <<
    //        num_flow_before_sparser << std::endl << std::endl;

    py_make_flow_sparser(num_in*num_out, min_flow);

    //int num_flow_after_sparser = get_num_nonzero_flow(min_flow);
    //std::cout << "after py_make_flow_sparser, num non-zero flow " <<
    //        num_flow_after_sparser << std::endl;


#ifdef LOG_SF
    shc_log_info(shc_logname, "after path decomm, num_flow %d\n", min_num_flow);
#endif

    if(flow_times > 1)
        return false;
    else
        return true;
}

void Sparse_flow_handler::
py_make_flow_sparser(int product, std::vector<double> & flow)
{
    if(product > sf_setting.path_sparsity)
    {
        std::vector<index_flow_t> index_flow;
        for(int i=0; i<flow.size(); i++)
            index_flow.emplace_back(i, flow[i]);

        Index_flow_sorter if_sorter;
        std::sort(index_flow.begin(), index_flow.end(), if_sorter);

        //for(int i=0; i<index_flow.size(); i++)
        //{
        //    modify_writer << index_flow[i].first << "\t" << index_flow[i].second << std::endl;
        //}


        flow.assign(flow.size(), 0);
        for(int i=0; i<sf_setting.path_sparsity; i++)
        {
            flow[index_flow[i].first] = index_flow[i].second;
        }
    }
}

void Sparse_flow_handler::remove_decomposed_node(std::set<vd_t> & remove_vd)
{
    for(std::set<vd_t>::iterator it  = remove_vd.begin();
                                    it != remove_vd.end(); it++)
    {
        boost::remove_vertex(*it, graph);
    }
}

bool is_connect_start_node(vd_t vd)
{

}
bool is_connect_end_node(vd_t vd)
{

}

void Sparse_flow_handler::update_sorted_vds(bool is_reverse)
{
    std::vector<vd_t> new_sorted_vds;
    for(std::vector<vd_t>::iterator it = sorted_vds.begin();
                                    it!= sorted_vds.end(); it++)
    {
        vd_t vd = *it;
        int in_degree = boost::in_degree(vd, graph);
        int out_degree = boost::out_degree(vd, graph);
        //shc_log_info(shc_logname, "node %u, in:%d out:%d\n", )
        if( (in_degree!=0 && out_degree!=0) &&
            vd!=vd_start && vd!=vd_end )
        {
            if(in_degree>1 || out_degree>1)
            {
                if (in_degree<=1)
                    continue;
                new_sorted_vds.push_back(vd);
                //shc_log_info(shc_logname, "vd %d is kept\n", graph[vd].node_id);
                log_local_node(vd, false);
            }
        }
    }
    //std::cout << new_sorted_vds.size() << "nodes are kept\n";
    //shc_log_info(shc_logname, "%d nodes are kept\n", new_sorted_vds.size());

    sorted_vds = new_sorted_vds;
    if(is_reverse)
        std::reverse(sorted_vds.begin(), sorted_vds.end());
}

void Sparse_flow_handler::modify_local_graph(vd_t v, std::vector<double> & flow,
                std::vector<vd_t> & in_nodes, std::vector<vd_t> & out_nodes)
{
    //shc_log_info(shc_logname, "Start modify graph\n");

    int num_in = boost::in_degree(v, graph);
    int num_out = boost::out_degree(v, graph);
    bundled_node_p & node_curr = graph[v];

    //modify_writer << "modify vd: " << node_curr.node_id
    //              << " in " << num_in << " out " << num_out << std::endl;
    //for(int k=0; k<flow.size(); k++)
    //{
    //    int i, j;
    //    get_i_j(k, num_out, i, j);
        //<< std::setprecision(2)
    //    modify_writer << "f(" << i << "," << j << ") = " << flow.at(k) << std::endl;
    //}

    for(int k=0; k<flow.size(); ++k)
    {
        int i, j;
        get_i_j(k, num_out, i, j);
        // use notation from paper
        if(flow[k] != 0)
        {
            // collecting information around v
            vd_t u = in_nodes[i];
            bundled_node_p & node_u = graph[u];
            aer_t u_v_aer = boost::edge(u, v, graph);
            bundled_edge_p & edge_u_v = graph[u_v_aer.first];

            vd_t w = out_nodes[j];
            bundled_node_p & node_w = graph[w];
            aer_t v_w_aer = boost::edge(v, w, graph);
            bundled_edge_p & edge_v_w = graph[v_w_aer.first];

            // create new nodes
            vd_t v_new = boost::add_vertex(graph);
            graph[v_new] = bundled_node_p(node_curr.seq, glob_node_id++);
            //std::cout << "graph[v_new].id " << graph[v_new].node_id << std::endl;
            //new_to_ori[v_new] = v;

            vd_contains_map.insert(std::make_pair(v_new, vd_contains_map[v]));

            aer_t aer_u_vn = boost::add_edge(u, v_new, graph);
            graph[aer_u_vn.first] = bundled_edge_p(edge_u_v.weight, flow[k]);
            aer_t aer_vn_w = boost::add_edge(v_new, w, graph);
            graph[aer_vn_w.first] = bundled_edge_p(edge_v_w.weight, flow[k]);
        }
    }
    boost::clear_vertex(v, graph);
    //shc_log_info(shc_logname, "Finish modify graph\n");
}

/** link nodes to front nodes and end nodes. exclude terminal nodes
 */
void Sparse_flow_handler::link_node_to_terminals(vd_t vd)
{
    if (vd==vd_end || vd==vd_start)
        return;

    if(boost::out_degree(vd, graph)==0)
    {
        aer_t aer = boost::add_edge(vd, vd_end, graph);
        graph[aer.first] = bundled_edge_p(0, 0);
    }

    if(boost::in_degree(vd, graph)==0)
    {
        //shc_log_info(shc_logname, "add leaf edge\n");
        aer_t aer = boost::add_edge(vd_start, vd, graph);
        graph[aer.first] = bundled_edge_p(0, 0);
    }
}

void Sparse_flow_handler::link_leaf_node_to_terminal(
                std::vector<vd_t> & in_nodes, std::vector<vd_t> & out_nodes)
{
    aer_t aer;
    for(std::vector<vd_t>::iterator it = in_nodes.begin();
                                    it!= in_nodes.end(); it++)
    {
        if(boost::out_degree(*it, graph)==0)
        {
            //shc_log_info(shc_logname, "add leaf edge\n");
            aer = boost::add_edge(*it, vd_end, graph);
            graph[aer.first] = bundled_edge_p(0, 0);
        }
    }

    for(std::vector<vd_t>::iterator it = out_nodes.begin();
                                    it!= out_nodes.end(); it++)
    {
        if(boost::in_degree(*it, graph)==0)
        {
            //shc_log_info(shc_logname, "add leaf edge\n");
            aer = boost::add_edge(vd_start, *it, graph);
            graph[aer.first] = bundled_edge_p(0, 0);
        }
    }
}

void Sparse_flow_handler::
make_in_out_flow_equal(std::vector<double> & in_counts,
                       std::vector<double> & out_counts)
{
    double in_sum = std::accumulate(in_counts.begin(), in_counts.end(), 0.0);
    double out_sum = std::accumulate(out_counts.begin(), out_counts.end(), 0.0);
    if(in_sum > out_sum)
    {
        double diff = in_sum - out_sum;
        std::for_each(out_counts.begin(), out_counts.end(),
             [diff, out_sum](double & k)
             {
                 k = (1.0+diff/out_sum) * k;
             }
        );
    }
    else
    {
        double diff = out_sum - in_sum;
        std::for_each(in_counts.begin(), in_counts.end(),
             [diff, in_sum](double & k)
             {
                 k = (1.0+diff/in_sum) * k;
             }
        );
    }
}

bool Sparse_flow_handler::
is_update_flow(std::vector<double> & flow, std::vector<double> & min_flow,
                std::vector<int>& sub_flow_indicator, double min_unknown_flow)
{
    bool cond_1 = min_flow.size()==0;
    //bool cond_2 = get_flow_difference(flow, min_flow) < sf_setting.tol &&
    //              dot_product(sub_flow_indicator, flow) < min_unknown_flow;
    double flow_diff = get_flow_difference(flow, min_flow);
    //double flow_dot_product = dot_product(sub_flow_indicator, flow);
    //modify_writer << "flow_diff " << flow_diff << " flow_dot_product " << flow_dot_product << std::endl;
    //modify_writer << "sf_setting.tol " << sf_setting.tol << "min_unknown_flow: " << min_unknown_flow << std::endl;
    bool cond_2 =  flow_diff < sf_setting.tol;// &&
                //flow_dot_product < min_unknown_flow;
    //bool cond_3 = std::accumulate(flow.begin(), flow.end(), 0.0) >
    //              std::accumulate(min_flow.begin(), min_flow.end(), 0.0);

    //double a =    std::accumulate(flow.begin(), flow.end(), 0.0);
    //double b =    std::accumulate(min_flow.begin(), min_flow.end(), 0.0);
    /*
    modify_writer << "a " << a << " b " << b << std::endl;
    if(cond_1 )
        modify_writer << "cond1" << std::endl;
    if(cond_2 )
        modify_writer << "cond2" << std::endl;
    if(cond_3 )
        modify_writer << "cond3" << std::endl;
    */
    return cond_1 || cond_2;
}

/**
 * This function truncate flow solution after LP, flow and flow_2 are copy of
 * each other, m is num of in edges, n is for num of out edges,
 */
void Sparse_flow_handler::truncate_flow_value(int m, int n,
                  std::vector<double> & flow, std::vector<double> & flow_2,
        std::vector<double> & in_counts, std::vector<double> & out_counts)
{
    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)
        {
            double sub_flow = flow[i*n+j];
            double min_count = std::min(in_counts[i], out_counts[j]);
            if (sub_flow < sf_setting.removal_factor*min_count)
            {
                bool is_only_out_flow = true;
                bool is_out_appear_once = false;
                // check if other flow passing this in node
                for(int k=0; k<n; k++)
                {
                    if(flow[i*n+k] > 0)
                    {
                        if(is_out_appear_once) // flow out
                        {
                            is_only_out_flow = false;
                            break;
                        }
                        is_out_appear_once = true;
                    }
                }
                if(is_only_out_flow)
                    continue;
                bool is_only_in_flow = true;
                bool is_in_appear_once = false;
                for(int k=0; k<m; k++)
                {
                    if(flow[k*n+j] > 0)
                    {
                        if(is_in_appear_once)
                        {
                            is_only_in_flow = false;
                            break;
                        }
                        is_in_appear_once = true;
                    }
                }
                if(is_only_in_flow)
                    continue;

                flow[i*n+j] = 0;
            }
        }
    }
}

void Sparse_flow_handler::old_truncate_flow_value(int m, int n,
             std::vector<double> & flow, std::vector<double> & flow_2,
             std::vector<double> & in_counts, std::vector<double> & out_counts)
{
    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)
        {
            double min_count = std::min(in_counts[i], out_counts[j]);
#ifdef LOG_LP_SUMMARY
            shc_log_info(shc_logname, "in %d: %f, out %d: %f\n", i, in_counts[i], j, out_counts[j]);
            shc_log_info(shc_logname, "flow %f tol %f rm %f\n",
                flow[i*n+j], sf_setting.tol, sf_setting.removal_factor*min_count);
#endif

            // for flow 2
            if (flow_2[i*n+j] < sf_setting.sparsity_factor*min_count)
            {
                flow_2[i*n+j] = 0;
            }

            if (flow[i*n+j] < sf_setting.removal_factor*min_count ||
                flow[i*n+j] < sf_setting.tol)
            {
#ifdef LOG_LP_SUMMARY
                shc_log_info(shc_logname, "flow %d %d set to 0\n", i, j);
#endif
                flow[i*n+j] = 0;
                flow_2[i*n+j] = 0;
            }

            if (flow_2[i*n+j] < 0)
            {
                flow[i*n+j] = 0;
                flow_2[i*n+j] = 0;
            }
#ifdef LOG_LP_SUMMARY
            shc_log_info(shc_logname, "\n");
#endif
            if(flow[i*n+j] < 0)
            {
                shc_log_error("LP optimize, get negative flow\n");
                exit(1);
            }
        }
    }
}

double Sparse_flow_handler::path_decompose(vd_t vd, std::vector<double>& gamma,
        std::vector<double> & equ_constraint,
        std::vector<int>& sub_flow_indicator, std::vector<double> & flow, double & scale)
{
    //shc_log_info(shc_logname, "Start decompose\n");
    glp_prob *lp;
    glp_term_out(GLP_OFF);
    std::vector<int> ia(1,0);
    std::vector<int> ja(1,0);
    std::vector<double> ar(1,0.0);

    double obj_value;

    int num_in = boost::in_degree(vd, graph);
    int num_out = boost::out_degree(vd, graph);

    int num_cost_var = num_in * num_out;
    int num_equ_constraint = num_in + num_out;

    assert(gamma.size()==num_cost_var);

    //create problem
    lp = glp_create_prob();
    glp_set_prob_name(lp, "SF LP");
    glp_set_obj_dir(lp, GLP_MIN);

    // fill problem
    glp_add_rows(lp, num_equ_constraint);
    for(int i=1; i<=num_equ_constraint; i++)
    {
        glp_set_row_bnds(lp, i, GLP_FX, equ_constraint.at(i-1),
                                        equ_constraint.at(i-1));
    }

    // construct objective
    glp_add_cols(lp, num_cost_var);
    for(int i=1; i<=num_cost_var; i++)
    {
        glp_set_col_bnds(lp, i, GLP_LO, 0, IGNORE); //f(e,e) greater than 0
        glp_set_obj_coef(lp, i, gamma.at(i-1)*sub_flow_indicator.at(i-1));
    }

    // construct constraint coeff
    int num_coeff = 0;
    for(int i=1; i<=num_in; i++)
    {
        for(int j=1; j<=num_out; j++)
        {
            ia.push_back(i);
            ja.push_back(num_out*(i-1) + j);
            ar.push_back(1.0);
            num_coeff++;
        }
    }
    // for out constraint
    for(int i=1; i<=num_out; i++)
    {
        for(int j=1; j<=num_in ; j++)
        {
            ia.push_back(i + num_in);
            ja.push_back(i + (j-1)*num_out);
            ar.push_back(1.0);
            num_coeff++;
        }
    }
#ifdef LOG_LP_SUMMARY
    shc_log_info(shc_logname, "\t\tproblem summary on vd %d\n", graph[vd].node_id);
    shc_log_info(shc_logname, "\t\tinD: %d, outD: %d, num unknown %d\n",
                            num_in, num_out, sum_vector(sub_flow_indicator));

    shc_log_info(shc_logname, "%d OBJ COST\n", gamma.size());
    for(int i=1; i<=num_cost_var; i++)
    {
        shc_log_info(shc_logname, "%d: %f\n", i,
                                    gamma.at(i-1)*sub_flow_indicator.at(i-1));
    }
    shc_log_info(shc_logname, "\t\t%d Equal constraint\n", equ_constraint.size());
    for(int k=1; k<=num_equ_constraint; k++)
    {
        shc_log_info(shc_logname, "%d: xxx = %f\n", k, scale*equ_constraint.at(k-1));

    }
    for(int k=1; k<=num_coeff; k++)
    {
        shc_log_info(shc_logname, "ia[%d] = %d, ja[%d] = %d, ar[%d] = %f\n",
                                k,ia.at(k), k,ja.at(k), k,ar.at(k));
    }
#endif

    // solve
    assert(ia.size()==ja.size());
    glp_load_matrix(lp, num_coeff, &ia.at(0), &ja.at(0), &ar.at(0));
    glp_simplex(lp, NULL);
    obj_value = glp_get_obj_val(lp);

    for(int i=1; i<=num_cost_var; i++)
        flow.push_back(glp_get_col_prim(lp, i));

    for(std::vector<double>::iterator it=flow.begin();
                           it!=flow.end(); it++)
    {
       (*it) *= scale;
       //std::cout << *it << std::endl;
    }

#ifdef LOG_LP_SUMMARY
    shc_log_info(shc_logname, "solution summary\n");
    for(int k=0; k<flow.size(); k++)
    {
        int i, j;
        get_i_j(k, num_out, i, j);
        shc_log_info(shc_logname, "f(%d,%d) = %f\n", i,j, flow.at(k));
    }
#endif
    glp_delete_prob(lp);
    glp_free_env();

    //shc_log_info(shc_logname, "Finish decompose\n");
    return obj_value;
}


void Sparse_flow_handler::
dump_SF_seq(std::string & out_file)
{
#ifdef LOG_SF
    shc_log_info(shc_logname, "start reconstruct seq\n");
#endif
    std::string curr_seq;
    std::vector<int> delimit;
    std::vector<vd_t> stack_vd;
    stack_vd.push_back(vd_start);

    std::ofstream file_writer(out_file);
    std::cout << "vd end in degree " << boost::in_degree(vd_end, graph)  << std::endl;
    std::cout << "vd start out degree " << boost::out_degree(vd_start, graph)  << std::endl;
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd=*it;
        if(vd != vd_end)
        {
            if(boost::in_degree(vd, graph) > 1)
            {
                shc_log_error("node %u has indegree > 1\n", graph[vd].node_id);
                log_node_info(vd, true);
                exit(0);
            }
        }
    }

    log_all_node(true);

    output_seq_helper(file_writer, curr_seq, delimit, stack_vd);

    file_writer.close();
}

void Sparse_flow_handler::
dump_SF_seq(std::string & out_file, pthread_mutex_t * writer_lock_ptr)
{
#ifdef LOG_SF
    shc_log_info(shc_logname, "start reconstruct seq multi thread\n");
#endif
    std::string curr_seq;
    std::vector<int> delimit;
    std::vector<vd_t> stack_vd;
    stack_vd.push_back(vd_start);

    std::ofstream file_writer;
    file_writer.open(out_file, std::ofstream::out | std::ofstream::app);

    output_seq_helper_mt(file_writer, writer_lock_ptr, curr_seq, delimit, stack_vd);

    file_writer.close();
}

void Sparse_flow_handler::
output_seq_helper(std::ofstream & file_writer, std::string & curr_seq,
                std::vector<int> & delimit, std::vector<vd_t> & stack_vd)
{
    vd_t vd = stack_vd.back();
    //std::cout << "VD: " << graph[vd].node_id << std::endl;
    //shc_log_info(shc_logname, "VD: %u\n", graph[vd].node_id);
    if( vd == vd_end)
    {
        std::string out_seq = curr_seq.substr(0, delimit.at(delimit.size()-2));
        if(out_seq.size() > setting.output_seq_min_len)
        {
            file_writer << ">comp_" << comp_id << "_graph_" << graph_id << "_"
                        << (num_out_seq++) << std::endl;
            file_writer << out_seq << std::endl;
        }
    }
    else
    {
        out_eip_t out_eip = boost::out_edges(vd, graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            bundled_edge_p & edge_p = graph[*out_ei];
            curr_seq += graph[vd_target].seq.substr(edge_p.weight);
            delimit.push_back(curr_seq.size());
            stack_vd.push_back(vd_target);

            output_seq_helper(file_writer, curr_seq, delimit, stack_vd);

            stack_vd.pop_back();
            delimit.pop_back();
            curr_seq.resize(delimit.back());
        }
    }
}

void Sparse_flow_handler::
output_seq_helper_mt(std::ofstream & file_writer, pthread_mutex_t * writer_lock_ptr,
        std::string & curr_seq, std::vector<int> & delimit, std::vector<vd_t> & stack_vd)
{
    vd_t vd = stack_vd.back();
    //std::cout << "VD: " << graph[vd].node_id << std::endl;
    //shc_log_info(shc_logname, "VD: %u\n", graph[vd].node_id);
    if( vd == vd_end)
    {
        std::string out_seq = curr_seq.substr(0, delimit.at(delimit.size()-2));
        if(out_seq.size() > setting.output_seq_min_len)
        {
            pthread_mutex_lock(writer_lock_ptr);
            file_writer << ">comp_" << comp_id << "_graph_" << graph_id << "_"
                        << (num_out_seq++) << std::endl;
            file_writer << (curr_seq.substr(0, delimit.at(delimit.size()-2))) << std::endl;
            pthread_mutex_unlock(writer_lock_ptr);
        }
    }
    else
    {
        if(boost::out_degree(vd, graph) > 0)
        {
            out_eip_t out_eip = boost::out_edges(vd, graph);
            for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
            {
                vd_t vd_target = boost::target(*out_ei, graph);
                bundled_edge_p & edge_p = graph[*out_ei];
                curr_seq += graph[vd_target].seq.substr(edge_p.weight);
                delimit.push_back(curr_seq.size());
                stack_vd.push_back(vd_target);

                output_seq_helper_mt(file_writer, writer_lock_ptr, curr_seq, delimit, stack_vd);

                stack_vd.pop_back();
                delimit.pop_back();
                curr_seq.resize(delimit.back());
            }
        }
    }
}

/**
 * This function constructs known path indicator which is a m*n vectors,
 * indicating which f(e,e' ) is absolutely in the flow, by setting 0 for that
 * index. The format is as following, for 3 in edges, 2 out edges
 * [f(1,1), f(1,2), f(2,1), f(2,2), f(3,1), f(3,2)]
 */
void Sparse_flow_handler::
get_sub_flow_indicator(vd_t vd, std::vector<int> & sub_flow_indicator,
                    std::vector<vd_t> & in_nodes, std::vector<vd_t> & out_nodes)
{
#ifdef LOG_SF
    shc_log_info(shc_logname, "Start sub flow indicator\n");
#endif
    // get all paths for this node
    Vd_Vds_map_iterator vd_conatins_it = vd_contains_map.find(vd);
    if(vd_conatins_it == vd_contains_map.end())
    {
        shc_log_warning("contain map does not conatin\n");
        exit(0);
    }
    std::deque<vd_t> & contained_vds = vd_conatins_it->second;

    for(int i=0; i<in_nodes.size(); i++)
    {
        vd_t vd_source = in_nodes[i];
        vd_t u = vd_contains_map[vd_source].back();
        for(int j=0; j<out_nodes.size(); j++)
        {
            vd_t vd_target = out_nodes[j];
            vd_t w = vd_contains_map[vd_target].front();

            if(!is_subflow_on_vd_unknown(contained_vds.front(), u, w, 1, contained_vds.size()))
            {
                sub_flow_indicator.push_back(0);
                continue;
            }

            if(!is_subflow_on_vd_unknown(contained_vds.back(), u, w, contained_vds.size(), 1))
            {
                sub_flow_indicator.push_back(0);
                continue;
            }

            // push the known or unknown to indicator
            sub_flow_indicator.push_back(1);
        }
    }

#ifdef LOG_SF
    for(int i=0; i< sub_flow_indicator.size(); i++)
    {
        shc_log_info(shc_logname, "%d flow ind is %d\n", i, sub_flow_indicator[i]);
    }
#endif
}

/**
 * check two conditions around vd, which is an original node loaded from file,
 * for a single path which contains vd.
 * 1. if u is b_dist before location of vd inside a path ,
 * 2. if w is a_dist after location of vd inside a path ,
 * return false, if both conditions are true, i.e. this vd is not unknown,hence known
 * return true otherwise
 */
bool Sparse_flow_handler::
is_subflow_on_vd_unknown(vd_t vd, vd_t u, vd_t w, int b_dist, int a_dist)
{
#ifdef LOG_SF
    shc_log_info(shc_logname, "vd:%d, L:%d, R:%d, LD:%d, RD:%d\n", graph[vd].node_id,
                graph[u].node_id, graph[w].node_id, b_dist, a_dist);
#endif
    Node_path_MMap_iterator_pair path_it_pair =
                                    paths_for_node.equal_range(vd);
    for(Node_path_MMap_iterator path_it = path_it_pair.first;
                                path_it!= path_it_pair.second; path_it++)
    {
        Path_info & path_info = path_it->second;

        //shc_log_info(shc_logname, "Path %d, location %d\n", path_info.id, path_info.location);

        std::vector<vd_t> & path = known_paths[path_info.id];
        if ( path_info.location-b_dist<0 ||
             path_info.location+a_dist >=path.size()-1 )
        {
            continue; // check next path
        }

        vd_t vd_before = path[path_info.location-b_dist];
        vd_t vd_after = path[path_info.location+a_dist];

        if(vd_before==u && vd_after==w)
        {
#ifdef LOG_SF
            shc_log_info(shc_logname, "Find vd, Path %d, location %d\n", path_info.id, path_info.location);
#endif
            return false;
        }
    }
#ifdef LOG_SF
    shc_log_info(shc_logname, "Find No\n");
#endif
    return true;
}
/*
vd_t Sparse_flow_handler:: get_ori_vd(vd_t vd)
{
    Vd_Vd_map_iterator it = new_to_ori.find(vd);
    if(it == new_to_ori.end())
        return vd;
    else
        return it->second;
}
*/
/**
 * Used for checking final answer if the entire graph is acyclic
 * @return false when not acyclic
 */
bool Sparse_flow_handler::acyclic_check()
{
    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        S_info & dfs_info = graph[vd].s_info;
        dfs_info.color = WHITE;
        dfs_info.parent = NULL;
    }

    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        S_info & dfs_info = graph[vd].s_info;
        if(dfs_info.color == WHITE)
        {
            if(!DFS_visit(vd))
            {
                shc_log_info(shc_logname, "cycle node vd %d\n", graph[vd].node_id);
                return false;
            }
        }
    }
    return true;
}

// return false when it is not acyclic
// implement from <<introduction to algorithm>>
bool Sparse_flow_handler::DFS_visit(vd_t vd)
{
    S_info & dfs_info = graph[vd].s_info;
    dfs_info.color = GRAY;
    out_eip_t out_eip = boost::out_edges(vd, graph);
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        S_info & v_dfs_info = graph[vd_target].s_info;
        if(v_dfs_info.color == GRAY)
        {
            shc_log_info(shc_logname, "cycle node vd %d\n", graph[vd].node_id);
            return false;
        }

        if(v_dfs_info.color == WHITE)
        {
            v_dfs_info.parent = vd;
            if(!DFS_visit(vd_target))
            {
                shc_log_info(shc_logname, "cycle node vd %d\n", graph[vd].node_id);
                return false;
            }
        }
    }
    dfs_info.color = BLACK;
    return true;
}

/**
 * Simple helper function, which puts in, out edges weight into vectors
 */
void Sparse_flow_handler::
get_in_out_edges(vd_t vd, std::vector<vd_t> & in_nodes, std::vector<vd_t> & out_nodes,
                 std::vector<double> & in_counts, std::vector<double> & out_counts)
{
    bundled_node_p & node = graph[vd];
    //shc_log_info(shc_logname, "node %u, inD: %d, outD: %d\n", node.node_id,
    //                boost::in_degree(vd,graph), boost::out_degree(vd,graph));
    in_eip_t in_eip = boost::in_edges(vd,graph);
    for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
    {
        vd_t vd_source = boost::source(*in_ei, graph);
        in_nodes.push_back(vd_source);
        bundled_edge_p & edge_p = graph[*in_ei];
        in_counts.push_back(edge_p.count);
    }
    out_eip_t out_eip = boost::out_edges(vd,graph);
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        out_nodes.push_back(vd_target);
        bundled_edge_p & edge_p = graph[*out_ei];
        out_counts.push_back(edge_p.count);
    }
#ifdef LOG_SF
    shc_log_info(shc_logname, "Finish get in out edge\n");
#endif
}

double Sparse_flow_handler::get_norm_2(std::vector<double> & v)
{
    return std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
}

double Sparse_flow_handler::get_norm_1(std::vector<double> & v)
{
    return std::accumulate(v.begin(), v.end(), 0.0, norm_1_helper);
}

int Sparse_flow_handler::get_norm_0(std::vector<double> & v)
{
    int count = 0;
    for(int i=0; i<v.size(); i++)
    {
        if(v[i] != 0)
            count++;
    }
    return count;
}

int Sparse_flow_handler::get_num_nonzero_flow(std::vector<double> & flow,
                                    std::vector<int> & sub_flow_indicator)
{
    int num = 0;
    assert(flow.size()==sub_flow_indicator.size());
    for(int i=0; i<sub_flow_indicator.size(); i++)
    {
        if(sub_flow_indicator[i] > 0 && flow[i]!=0)
        {
            num++;
        }
    }
    return num;
}

int Sparse_flow_handler::get_num_nonzero_flow(std::vector<double> & flow)
{
    int num = 0;
    for(int i=0; i<flow.size(); i++)
    {
        if(flow[i]!=0)
        {
            num++;
        }
    }
    return num;
}

double Sparse_flow_handler::get_flow_difference(std::vector<double> & flow1,
                               std::vector<double> & flow2)
{
    assert(flow1.size() == flow2.size());
    double sum_1 = std::accumulate(flow1.begin(), flow1.end(), 0.0);
    double sum_2 = std::accumulate(flow2.begin(), flow2.end(), 0.0);
    return abs(sum_1 - sum_2);
}

void Sparse_flow_handler::get_i_j(int flow_index, int num_out, int & i, int & j)
{
    i = flow_index/num_out;
    j = flow_index%num_out;
}

void Sparse_flow_handler::set_trials(int m, int n)
{
    sf_setting.num_trials = sf_setting.multiple_test*std::min(2*m*n*std::max(m,n), 100);
}

double Sparse_flow_handler::sum_vector(std::vector<double> & v)
{
    return std::accumulate(v.begin(), v.end(), 0.0);
}

int Sparse_flow_handler::sum_vector(std::vector<int> & v)
{
    return std::accumulate(v.begin(), v.end(), 0);
}

/**
 * for nodes like
 *  x \                         /  x
 *  x -  x - x but not[ x -  x -  x ]
 *  x /                        \  x
 *  merge middle nodes to the left nodes or right nodes respectively
 * return if delete this node
 */
bool Sparse_flow_handler::link_Y_path(vd_t vd)
{
    //shc_log_info(shc_logname, "\t\t link Y path\n");
    int in_degree = boost::in_degree(vd, graph);
    int out_degree = boost::out_degree(vd, graph);
    assert((out_degree==1 && in_degree>1));

    bundled_node_p & curr_node = graph[vd];
    // link in to out

    out_eip_t out_eip = boost::out_edges(vd,graph);
    vd_t vd_target = boost::target(*out_eip.first, graph);
    bundled_edge_p & out_edge_p = graph[*out_eip.first];

    //if(vd_target == vd_end)
    //    return false;

    double sum_in_count = 0;
    in_eip_t in_eip = boost::in_edges(vd,graph);
    for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
    {
        bundled_edge_p & edge_p = graph[*in_ei];
        sum_in_count += edge_p.count;
    }

    double scale = out_edge_p.count/sum_in_count;

    in_eip = boost::in_edges(vd,graph);
    for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
    {
        vd_t vd_source = boost::source(*in_ei, graph);
        bundled_node_p & source_node = graph[vd_source];
        bundled_edge_p & edge_p = graph[*in_ei];
        source_node.seq = source_node.seq + curr_node.seq.substr(edge_p.weight);
        aer_t aer = boost::add_edge(vd_source, vd_target, graph);
        graph[aer.first] = bundled_edge_p(out_edge_p.weight, edge_p.count*scale); //out_edge_p;

        //for path indicator
        update_contain_map(vd_source, vd);
    }

    //clean vertex
    //shc_log_info(shc_logname, "linked Y path, node %u, inD %d, outD %d\n",
    //                            curr_node.node_id, in_degree, out_degree);
    boost::clear_vertex(vd, graph);
    //shc_log_info(shc_logname, "after clear vertex\n");
    return true;
}

void Sparse_flow_handler::update_contain_map(vd_t vd_containing, vd_t vd_contained)
{
    //for path indicator
    Vd_Vds_map_iterator vd_conatins_it = vd_contains_map.find(vd_containing);
    if(vd_conatins_it !=  vd_contains_map.end())
    {
        Vd_Vds_map_iterator vd_contained_it = vd_contains_map.find(vd_contained);
        if(vd_contained_it==vd_contains_map.end())
        {
            shc_log_warning("contain does not have vd\n");
            exit(0);
        }
        else
        {
            std::deque<vd_t> & vec_of_contained = vd_contained_it->second;
            std::deque<vd_t>::iterator it_of_containing = vd_conatins_it->second.begin();
            vd_conatins_it->second.insert(it_of_containing,
                            vec_of_contained.begin(), vec_of_contained.end());
        }
    }
    else
    {
        shc_log_warning("contain does not have vd\n");
        exit(0);
    }
}

// debug
void Sparse_flow_handler::log_graph_struct (bool show_seq){
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        shc_log_info(shc_logname, "******************************************\n");
        if (vd == vd_start)
            shc_log_info(shc_logname, "START NODE\n");
        if (vd == vd_end)
            shc_log_info(shc_logname, "END NODE\n");
        shc_log_info(shc_logname, "VD: %u\n", graph[vd].node_id);
        if(show_seq)
            shc_log_info(shc_logname, "Seq: %s\n", graph[vd].seq.c_str());
        in_eip_t in_eip = boost::in_edges(vd,graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
            vd_t vd_source = boost::source(*in_ei, graph);
            shc_log_info(shc_logname, "IN  VD: %6u, weight %6d, count %4.2f\n",
                graph[vd_source].node_id, graph[*in_ei].weight, graph[*in_ei].count);
        }
        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            shc_log_info(shc_logname, "OUT VD: %6u, weight %6d, count %4.2f\n",
                graph[vd_target].node_id, graph[*out_ei].weight, graph[*out_ei].count);
        }
    }
}

void Sparse_flow_handler::log_term_node(bool is_front, bool is_back, bool detail)
{
    if(is_front)
    {
        shc_log_info(shc_logname, "Log start node, inD: %d, outD: %d\n",
           boost::in_degree(vd_start,graph), boost::out_degree(vd_start,graph));
        in_eip_t in_eip = boost::in_edges(vd_start, graph);
        if(detail)
        {
            for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
            {
                vd_t vd_source = boost::source(*in_ei, graph);
                shc_log_info(shc_logname, "IN  VD: %6u, weight %6d, count %4.2d\n",
                    graph[vd_source].node_id, graph[*in_ei].weight, graph[*in_ei].count);
            }
        }
    }
    if(is_back)
    {
        shc_log_info(shc_logname, "Log end node, inD: %d, outD: %d\n",
           boost::in_degree(vd_end,graph), boost::out_degree(vd_end,graph));
        out_eip_t out_eip = boost::out_edges(vd_end,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            if(detail)
            {
                vd_t vd_target = boost::target(*out_ei, graph);
                shc_log_info(shc_logname, "OUT VD: %6u, weight %6d, count %4.2f\n",
                    graph[vd_target].node_id, graph[*out_ei].weight, graph[*out_ei].count);
            }
        }
    }
}

void Sparse_flow_handler::log_local_node(vd_t vd, bool is_show_seq)
{
    shc_log_info(shc_logname, "\t\t\t VD:%6d\n", graph[vd].node_id);
    if (is_show_seq)
        shc_log_info(shc_logname, "SEQ: %s\n", graph[vd].seq.c_str());

    in_eip_t in_eip = boost::in_edges(vd,graph);
    for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
    {
        vd_t vd_source = boost::source(*in_ei, graph);
        shc_log_info(shc_logname, "IN  VD: %6d\n", graph[vd_source].node_id);
    }
    out_eip_t out_eip = boost::out_edges(vd,graph);
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        shc_log_info(shc_logname, "OUT VD: %6d\n",graph[vd_target].node_id);
    }
    Vd_Vds_map_iterator contain_it = vd_contains_map.find(vd);
    std::deque<vd_t> & vec = contain_it->second;
    for(int i=0; i<vec.size();i++)
    {
        shc_log_info(shc_logname, "contain VD: %d\n", graph[vec[i]].node_id);
    }
}

void Sparse_flow_handler::log_around_local_node(vd_t vd,
                std::vector<vd_t> & in_nodes, std::vector<vd_t> & out_nodes)
{
    log_local_node(vd, true);

    for(std::vector<vd_t>::iterator it = in_nodes.begin();
                                    it!= in_nodes.end(); ++it)
    {
        shc_log_info(shc_logname, "left : %6d\n", graph[*it].node_id);
        out_eip_t out_eip = boost::out_edges(*it, graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            shc_log_info(shc_logname, "%6d\n", graph[vd_target].node_id);
        }
    }

    for(std::vector<vd_t>::iterator it = out_nodes.begin();
                                    it!= out_nodes.end(); ++it)
    {
        shc_log_info(shc_logname, "right: %6d\n", graph[*it].node_id);


        in_eip_t in_eip = boost::in_edges(*it, graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
            vd_t vd_source = boost::source(*in_ei, graph);
            shc_log_info(shc_logname, "%6d\n", graph[vd_source].node_id);
        }
    }
}


// This function takes one argument, condense left or right if it can,
// otherwise return false
bool Sparse_flow_handler::
condense(vd_t vd, bool is_check_left , vd_t vd_new)
{
    //shc_log_info(shc_logname, "condense local\n");
    vd_t vd_target, vd_source;
    if(is_check_left)
    {
        vd_target = vd;
        if(boost::in_degree(vd_target,graph)!=1)
        {
            return false;
        }
        else
        {
            in_eip_t in_eip = boost::in_edges(vd_target, graph);
            vd_source = boost::source(*(in_eip.first), graph);
            if (vd_source== vd_start)
                return false;
            if(boost::out_degree(vd_source, graph)!=1)
                return false;
            //update info
            bundled_node_p & new_node = graph[vd_new];
            bundled_node_p & source_node = graph[vd_source];
            new_node.count = (new_node.count + source_node.count)/2;
            bundled_edge_p & edge_p = graph[*(in_eip.first)];

            std::string seq_back;
            seq_back.resize(source_node.seq.size());
            std::reverse_copy(source_node.seq.begin(), source_node.seq.end(), seq_back.begin());
            new_node.seq += seq_back.substr(edge_p.weight);
        }
    }
    else
    {
        vd_source = vd;
        if(boost::out_degree(vd_source, graph)!=1)
        {
            return false;
        }
        else
        {
            out_eip_t out_eip = boost::out_edges(vd_source, graph);
            vd_target = boost::target(*(out_eip.first), graph);
            if(boost::in_degree(vd_target, graph)!=1)
                return false;
            if (vd_target== vd_end)
                return false;
            bundled_node_p & new_node = graph[vd_new];
            bundled_node_p & target_node = graph[vd_target];
            new_node.count = (new_node.count + target_node.count)/2;
            // info is only updated in the left side
            bundled_edge_p & edge_p = graph[*(out_eip.first)];
            int seq_len_before_concat = new_node.seq_len();
            new_node.seq = new_node.seq + target_node.seq.substr(edge_p.weight);
        }
    }

    assert(boost::out_degree(vd_source, graph) == 1 &&
           boost::in_degree(vd_target, graph) == 1);

    return true;
}

void Sparse_flow_handler::
condense_left(vd_t vd_target, vd_t vd_new, std::set<vd_t> & vd_remove_set)
{
    vd_t vd_source;
    while (condense(vd_target, true, vd_new))
    {
        //get the source node and clear
        in_eip_t in_eip = boost::in_edges(vd_target, graph);
        vd_source = boost::source(*(in_eip.first), graph);

        // for path indicator
        update_contain_map(vd_new, vd_source);

        boost::clear_in_edges(vd_target, graph);
        vd_remove_set.insert(vd_source);

        //prepare for next iteration
        vd_target = vd_source;
    }
    vd_t vd_leftmost = vd_target;
    if(vd_leftmost != vd_new)
    {    //link left nodes to new nodes
        in_eip_t in_eip = boost::in_edges(vd_leftmost, graph);
        for(in_ei_t it = in_eip.first; it!=in_eip.second; it++)
        {
            aer_t aer = boost::add_edge(boost::source(*it, graph), vd_new, graph);
            //copy the edge count
            graph[aer.first] = bundled_edge_p(graph[*it].weight, graph[*it].count);
        }
        boost::clear_in_edges(vd_leftmost, graph);
    }
}


void Sparse_flow_handler::
condense_right(vd_t vd_source, vd_t vd_new ,std::set<vd_t> & vd_remove_set)
{
    vd_t vd_target;
    while (condense(vd_source, false, vd_new))
    {
        //get the source node and clear
        out_eip_t out_eip = boost::out_edges(vd_source, graph);
        vd_target = boost::target(*(out_eip.first), graph);

        // for path indicator
        update_contain_map(vd_new, vd_target);

        boost::clear_out_edges(vd_source, graph);
        vd_remove_set.insert(vd_target);

        //prepare for next iteration
        vd_source = vd_target;
    }

    vd_t vd_rightmost = vd_source; // the node which stands for the right end
    //link right node
    if(vd_rightmost != vd_new)
    {
        out_eip_t out_eip = boost::out_edges(vd_rightmost, graph);
        for(out_ei_t it = out_eip.first; it!=out_eip.second; it++)
        {
            aer_t aer = boost::add_edge(vd_new, boost::target(*it, graph), graph);
            //copy the edge
            graph[aer.first] = bundled_edge_p(graph[*it].weight, graph[*it].count);
        }
        boost::clear_out_edges(vd_rightmost, graph);
    }
}

bool Sparse_flow_handler::is_condensable(vd_t vd, bool is_left)
{
    vd_t vd_source, vd_target;
    if(is_left)
    {
        if (boost::in_degree(vd, graph) != 1)
            return false;
        vd_target = vd;
        in_eip_t in_eip = boost::in_edges(vd, graph);
        vd_source = boost::source(*(in_eip.first), graph);
        if(vd_source == vd_start)
            return false;
        return (boost::out_degree(vd_source, graph) == 1);
    }
    else
    {
        if (boost::out_degree(vd, graph) != 1)
            return false;
        vd_source = vd;
        out_eip_t out_eip = boost::out_edges(vd, graph);
        vd_target = boost::target(*(out_eip.first), graph);
        if (vd_target == vd_end)
            return false;
        return (boost::in_degree(vd_target, graph) == 1);
    }
}

// always remember to clear the vd_remove_set afterward
vd_t Sparse_flow_handler::
local_condense(vd_t vd, std::set<vd_t> & vd_remove_set)
{
    bool left_condensable = is_condensable(vd, true);
    bool right_condensable = is_condensable(vd, false);
    //if(left_condensable)
    //    shc_log_info(shc_logname, "condense left\n");
    //if(right_condensable)
    //    shc_log_info(shc_logname, "condense right\n");
    if( left_condensable || right_condensable)
    {
        std::reverse(graph[vd].seq.begin(), graph[vd].seq.end());
        condense_left(vd, vd, vd_remove_set);
        std::reverse(graph[vd].seq.begin(), graph[vd].seq.end());
        condense_right(vd, vd ,vd_remove_set);
        return vd;
    }
    else
    {
        return NULL;//shc_log_info(shc_logname, "condense nothing\n");
    }
}

void Sparse_flow_handler::condense_graph(std::set<vd_t> & vd_remove_set)
{
    out_ei_t out_ei, out_ei_end;
    vi_t vi, vi_end, next;
    boost::tie(vi, vi_end) = boost::vertices(graph);
    //walk through the graph list, and update if
    uint64_t num_vertices_removed = 0;
    for(next=vi; vi!=vi_end; vi=next)
    {
        ++next;
        vd_t vd_source = *vi;
        vd_t vd_new = local_condense(vd_source, vd_remove_set);
        if(vd_new != NULL) // something condensed
        {
            for(std::set<vd_t>::iterator it=vd_remove_set.begin();
                it!=vd_remove_set.end(); ++it)
            {
                assert(boost::in_degree(*it, graph)==0 &&
                        boost::out_degree(*it, graph)==0);
            }
        }
        else
            shc_log_info(shc_logname, "no condense\n");
    }
}

void Sparse_flow_handler::check_any_condensible_node()
{
    shc_log_info(shc_logname, "check_any_condensible_node\n");
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        bool left_condensable = is_condensable(vd, true);
        bool right_condensable = is_condensable(vd, false);
        if(left_condensable)
        {
            shc_log_info(shc_logname, "VD: %d, left condensible\n",graph[vd].node_id);
        }
        if(right_condensable)
        {
            shc_log_info(shc_logname, "VD: %d, right condensible\n",graph[vd].node_id);
        }
    }
    shc_log_info(shc_logname, "Finish check_any_condensible_node\n");
}

void Sparse_flow_handler::log_all_node(bool is_log_node)
{
    shc_log_info(shc_logname, "log all nodes, number of vertices is %d, number of edges is %d\n",
                boost::num_vertices(graph), boost::num_edges(graph));
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        log_node_info(*it, is_log_node);
    }
}

void Sparse_flow_handler::log_node_info(vd_t vd, bool is_log_nodes)
{
    bundled_node_p & curr_node = graph[vd];
    assert(!curr_node.seq.empty());
    //assert(boost::in_degree(vd,graph)<=4 && boost::out_degree(vd,graph)<=4);
    //shc_log_info(shc_logname, "********************************************************\n");
    shc_log_info(shc_logname,"\n");
    shc_log_info(shc_logname,"\t\tLOG_NODE\n");
    shc_log_info(shc_logname, "Node ID: %d\n", curr_node.node_id);
    shc_log_info(shc_logname, "Seq: %s\n", curr_node.seq.c_str());
    shc_log_info(shc_logname, "in degree %d, out degree %d\n",
            boost::in_degree(vd,graph), boost::out_degree(vd,graph));
    if(is_log_nodes)
    {
        shc_log_info(shc_logname, "In nodes:\n");
        in_eip_t in_eip = boost::in_edges(vd,graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
            vd_t vd_source = boost::source(*in_ei, graph);
            shc_log_info(shc_logname, "%d %d %d %d %f %s\n", boost::in_degree(vd_source, graph),
                    boost::out_degree(vd_source, graph), graph[vd_source].node_id,
                    graph[*in_ei].weight, graph[*in_ei].count ,graph[vd_source].seq.c_str());
        }
        shc_log_info(shc_logname, "Out nodes:\n");
        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            shc_log_info(shc_logname, "%d %d %d %d %f %s\n", boost::in_degree(vd_target, graph),
                    boost::out_degree(vd_target, graph) , graph[vd_target].node_id,
                    graph[*out_ei].weight, graph[*out_ei].count, graph[vd_target].seq.c_str());
        }
    }
}
