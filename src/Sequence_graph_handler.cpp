#include <boost/graph/compressed_sparse_row_graph.hpp>

#include "Sequence_graph_handler.h"

void * Sequence_graph_handler::run_seq_graph_handler_helper()
{
    //seq_graph_handler.load_kmer_and_build_graph(comp_i);
    Block_timer main_timer;
    start_timer(&main_timer);
    build_kmer_graph_from_reads();

    condense_graph();

    update_xnode_set();

    assign_read_to_xnode();
    //seq_graph_handler.log_all_node(true, false);

    bridge_all_xnodes();

    break_self_loops();

    break_all_cycles();
    //shc_log_info(shc_logname,"QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n");
    find_known_path(max_hop);
    std::cout << "Job " << component_i << " finish, ";
    stop_timer(&main_timer);
}

void Sequence_graph_handler::assign_read_to_xnode()
{
    shc_log_info(shc_logname, "Start assign reads to x node\n");
    start_timer(&timer);
    std::string temp, read_base, kmer_base;

    kmer_base.resize(kmer_length);

    //construct a map from start seq kmer to vertex
    typedef std::multimap<std::string, vd_t> Seq_node_MMap;
    typedef std::multimap<std::string, vd_t>::iterator
                                             Seq_node_MMap_iterator;
    typedef std::pair<Seq_node_MMap_iterator, Seq_node_MMap_iterator>
                                        Seq_node_MMap_iterator_pair;
    Seq_node_MMap start_seq_node_mmap;

    shc_log_info(shc_logname, "constructing x-node k start seq map\n");
    for(std::set<vd_t>::iterator it=xnode_set.begin(); it!=xnode_set.end(); it++)
    {
        vd_t vd = *it;
        std::string k_start(graph[vd].seq.substr(0,kmer_length));
        start_seq_node_mmap.insert(std::make_pair(k_start,vd));
    }

    //examine each subseq
    shc_log_info(shc_logname, "Reading lines\n");
    for(int j = 0; j < read_list.get_num_read(); j++)
    {
        Read_acc acc = read_list.get_read(j);
        for(int i=0; i<acc.len - kmer_length + 1; i++)
        {
            memcpy(&kmer_base.at(0), acc.read_ptr + i, kmer_length);
            Seq_node_MMap_iterator_pair sn_it_pair =
                    start_seq_node_mmap.equal_range(kmer_base);

            for(Seq_node_MMap_iterator sn_it =sn_it_pair.first;
                                       sn_it!=sn_it_pair.second; ++sn_it)
            {
                if(is_read_bridge_node(acc.read_ptr, acc.len, i, sn_it->second))
                {
                    node_link_read(sn_it->second, j, i);
                }
            }
        }
    }
    shc_log_info(shc_logname, "Finish assign reads to x node\n");

#ifdef SHOW_SEQ_GRAPH_PROGRESS
    stop_timer(&timer);
    std::cout << "load and match reads finish, ";
#endif
}

void Sequence_graph_handler::
node_link_read(vd_t vd, read_num_t read_id, read_length_t start)
{

    //shc_log_info(shc_logname, "in degree %d, out degree %d\n",
    //        boost::in_degree(vd,graph), boost::out_degree(vd,graph));
    graph[vd].reads_info.emplace_back(read_id, start);
}

// check if loaded read seq covers base seq in the node, return true if Yes
bool Sequence_graph_handler::
is_read_bridge_node(std::string & read_seq, read_length_t start, vd_t vd)
{
    if (start==0)
        return false;
    if (read_seq.size() <= start + graph[vd].seq_len())
        return false;
    return memcmp(&graph[vd].seq.at(0),
            &read_seq.at(start), graph[vd].seq_len())==0;
}

// check if loaded read seq covers base seq in the node, return true if Yes
bool Sequence_graph_handler::
is_read_bridge_node(const char * read_ptr, read_length_t len,
                                        read_length_t info_start, vd_t vd)
{
    if (info_start==0)
        return false;
    if (len <= info_start + graph[vd].seq_len())
        return false;
    return memcmp(&graph[vd].seq.at(0),
            read_ptr+info_start, graph[vd].seq_len())==0;
}

// it allows self loops
void Sequence_graph_handler::build_kmer_graph_from_reads()
{
    start_timer(&timer);
    //load all kmer first
    shc_log_info(shc_logname, "Start load read kmer and build graph\n");

    std::ifstream file_reader(kmer_path.c_str());
    std::string count_str, kmer_base, temp, read_base;
    uint64_t byte = 0;
    kmer_count_t count;

    // create nodes and
    while(  std::getline(file_reader, kmer_base, '\t') &&
            std::getline(file_reader, count_str)   )
    {
        vd_t vd = boost::add_vertex(graph);
        count = std::stoi(count_str);
        graph[vd] = bundled_node_p(kmer_base, count);
        kmer_node_map[kmer_base] = vd;
    }
    file_reader.close();

    // start to load reads
    file_reader.open(read_path);
    shc_log_info(shc_logname, "Reading lines\n");
    kmer_base.resize(kmer_length);
    while(  std::getline(file_reader, temp) &&
            std::getline(file_reader, read_base))
    {
        vd_t prev_vd = NULL;
        for(int i=0; i<read_base.size()-kmer_length+1; i++)
        {
            memcpy(&kmer_base.at(0), &read_base.at(i), kmer_length);
            Kmer_Node_map_iterator knc_iter = kmer_node_map.find(kmer_base);
            // If there is such node
            if(knc_iter != kmer_node_map.end())
            {
                if (prev_vd != NULL)
                {
                    //shc_log_info(shc_logname, "extending\n");
                    aer_t aer = boost::edge(prev_vd, knc_iter->second, graph);
                    if(!aer.second)
                    {
                        aer = boost::add_edge(prev_vd, knc_iter->second, graph);
                        graph[aer.first] = bundled_edge_p(kmer_length-1);
                    }
                    //else
                    //{
                    //    graph[aer.first].copy_count++;
                    //}
                }
                prev_vd = knc_iter->second;
            }
            else
            {
                //shc_log_info(shc_logname, "broken read\n");
                prev_vd = NULL;
            }
        }

        read_list.add_read(read_base);
        read_base.clear();
    }
    shc_log_info(shc_logname, "Finish load read kmer and build graph\n");

#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "building graph finish, ";
    stop_timer(&timer);
#endif
}
/**
 * This function creates a new node which capture vd1->vd2
 * @param vd1
 * @param vd2
 */

Sequence_graph_handler::vd_t Sequence_graph_handler::
merge_two_vertices(vd_t vd1, vd_t vd2)
{
    //shc_log_info(shc_logname, "merging two vertices\n");
    aer_t aer;
    //create merged node
    aer = edge(vd1, vd2, graph);

    bundled_node_p & node1 = graph[vd1];
    bundled_node_p & node2 = graph[vd2];

    //shc_log_info(shc_logname, "node 1 %d\n", node1.seq.size());
    //shc_log_info(shc_logname, "node 2 %d\n", node2.seq.size());
    std::string right_str = node2.seq.substr(graph[aer.first].weight);
    std::string merged_seq = node1.seq + right_str;
    //shc_log_info(shc_logname, "right sizse %d\n", right_str.size());
    //shc_log_info(shc_logname, "merged size %d\n", merged_seq.size());

    kmer_count_t merged_count = (node1.count + node2.count)/2;
    vd_t vd = boost::add_vertex(graph);
    graph[vd] = bundled_node_p(merged_seq, merged_count);
    merged_seq.clear();

    //merge reads, if the reads are not empty
    // something here
    if(!node1.reads_info.empty())
    {
        graph[vd].reads_info = node1.reads_info;
    }

    //link in edges
    in_eip_t vd1_in_eip = boost::in_edges(vd1, graph);
    for(in_ei_t it = vd1_in_eip.first; it!=vd1_in_eip.second; it++)
    {
        aer = boost::add_edge(boost::source(*it, graph), vd, graph);
        graph[aer.first] = bundled_edge_p(graph[*it].weight);
    }
    //link out edges
    out_eip_t vd2_out_eip = boost::out_edges(vd2, graph);
    for(out_ei_t it = vd2_out_eip.first; it!=vd2_out_eip.second; it++)
    {
        aer = boost::add_edge(vd, boost::target(*it, graph), graph);
        graph[aer.first] = bundled_edge_p(graph[*it].weight);
    }

    //shc_log_info(shc_logname, "Finish merge vertices\n");
    return vd;
}

void Sequence_graph_handler::condense_graph()
{
    start_timer(&timer);
    shc_log_info(shc_logname, "Start condense graph\n");
    out_ei_t out_ei, out_ei_end;
    vi_t vi, vi_end, next;
    boost::tie(vi, vi_end) = boost::vertices(graph);
    //walk through the graph list, and update if
    shc_log_info(shc_logname, "before condense, number of vertices is %d\n", boost::num_vertices(graph));
    std::vector<vd_t> vd_remove_list;
    int num_vertices_removed = 0;
    for(next=vi; vi!=vi_end; vi=next)
    {
        ++next;
        vd_t vd_source = *vi;
        if(boost::out_degree(vd_source, graph)==1)
        {
            boost::tie(out_ei, out_ei_end) = boost::out_edges(vd_source, graph);

            vd_t vd_target = boost::target(*out_ei, graph);
            if(boost::in_degree(vd_target, graph)==1)
            {
                merge_two_vertices(vd_source, vd_target);
                vd_remove_list.push_back(vd_source);
                vd_remove_list.push_back(vd_target);
                boost::clear_vertex(vd_source,graph);
                boost::clear_vertex(vd_target,graph);
                num_vertices_removed++;
            }
        }
    }
    for(std::vector<vd_t>::iterator it=vd_remove_list.begin(); it != vd_remove_list.end(); it++)
        boost::remove_vertex(*it,graph);
    shc_log_info(shc_logname, "after condense, number of vertices is %d, number vertices removed is %d\n", boost::num_vertices(graph), num_vertices_removed);
    shc_log_info(shc_logname, "Finish condense graph\n");
#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "condense graph and remove nodes finish ";
    stop_timer(&timer);
#endif
}

void Sequence_graph_handler::update_xnode_set()
{
    start_timer(&timer);
    shc_log_info(shc_logname, "Start getting all xnodes\n");
    if(xnode_set.size()!=0)
        xnode_set.clear();
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        if(is_xnode(*it))
        {
            xnode_set.insert(*it);
        }
    }
    shc_log_info(shc_logname, "Finish getting all xnodes\n");
#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "get " << xnode_set.size() << " xnode out of " << boost::num_vertices(graph) <<" finish ";
    stop_timer(&timer);
#endif
}

bool Sequence_graph_handler::is_xnode(vd_t vd)
{
    return boost::out_degree(vd, graph)>=2 && boost::in_degree(vd, graph) >=2;
}

bool Sequence_graph_handler::
is_read_bridge_left(Bdg_read_info & read_info, vd_t vd)
{
    //shc_log_info(shc_logname,"check read bridge left\n");
    Read_acc acc = read_list.get_read(read_info.read_id);
    char left_read_base = *(acc.read_ptr + read_info.start - 1);
    in_eip_t in_eip = boost::in_edges(vd, graph);
    for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
    {
        vd_t vd_source = boost::source(*in_ei, graph);
        bundled_node_p & source_node = graph[vd_source];
        bundled_edge_p & in_edge = graph[*in_ei];

        if(source_node.seq_len()==in_edge.weight)
        {
            shc_log_info("encounter same length,check read bridge left\n");
            return is_read_bridge_left(read_info, vd_source);
        }

        char left_edge_base = source_node.seq.at(source_node.seq_len()-in_edge.weight-1);
        if( left_read_base == left_edge_base)
        {
            return true;
        }
    }
    return false;
}
bool Sequence_graph_handler::
is_read_bridge_right(Bdg_read_info & read_info, vd_t vd)
{
    //shc_log_info(shc_logname,"check read bridge right\n");
    Read_acc acc = read_list.get_read(read_info.read_id);
    bundled_node_p & curr_node = graph[vd];
    char right_read_base = *(acc.read_ptr + read_info.start + curr_node.seq_len());
    out_eip_t out_eip = boost::out_edges(vd, graph);
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        bundled_node_p & target_node = graph[vd_target];
        bundled_edge_p & out_edge = graph[*out_ei];
        //if the edge contain the whole overlap
        if(out_edge.weight==target_node.seq_len())
        {
            shc_log_info("encounter same length,check read bridge right\n");
            Bdg_read_info start_update_info(read_info.read_id,
                    read_info.start+curr_node.seq_len()-target_node.seq_len());
            return is_read_bridge_right(start_update_info, vd_target);
        }

        char right_edge_base = target_node.seq.at(out_edge.weight);
        if( right_read_base == right_edge_base)
        {
            return true;
        }
    }
    return false;
}

void Sequence_graph_handler::refresh_bridging_reads(vd_t vd)
{
    //shc_log_info(shc_logname, "refresh_bridging_reads\n");
    bundled_node_p & curr_node = graph[vd];
    std::vector<Bdg_read_info> & reads_info = curr_node.reads_info;
    std::vector<Bdg_read_info> valid_reads_info;

    for(std::vector<Bdg_read_info>::iterator it = reads_info.begin();
                                             it!= reads_info.end(); it++)
    {
        Read_acc acc = read_list.get_read(it->read_id);

        if(!is_read_bridge_node(acc.read_ptr, acc.len, it->start, vd))
            continue;

        if(is_read_bridge_left(*it, vd) && is_read_bridge_right(*it,vd))
        {
            valid_reads_info.push_back(*it);
        }
    }
    curr_node.reads_info.assign(valid_reads_info.begin(), valid_reads_info.end());
    //shc_log_info(shc_logname, "Finish refresh_bridging_reads\n");
}

/**
 * Case when no read bridge xnode
 *              x             x
 *                \          /
 *                 x-x-x-x-x
 *                /          \
 *              x              x
 * and no read covers horizontal part
 * @param vd
 * @return
 */
bool Sequence_graph_handler::is_xnode_bridged(vd_t vd)
{
    //shc_log_info(shc_logname, "checking xnode\n");
    bundled_node_p & curr_node = graph[vd];
    std::vector<Bdg_read_info> & reads_info = curr_node.reads_info;
    std::vector<Bdg_read_info>::iterator ri_it;

    //if there is no reads for that xnode
    if(reads_info.size()==0)
        return false;

    std::set<char> in_base_set;
    std::set<char> out_base_set;

    // if node bridging read is not empty
    //shc_log_info(shc_logname, "constructing in out base set\n");
    for(ri_it=reads_info.begin(); ri_it!=reads_info.end(); ++ri_it )
    {
        Read_acc acc = read_list.get_read(ri_it->read_id);

        //shc_log_info(shc_logname, "Node %s\n", graph[vd].seq.c_str());
        //log_read_seq(acc);
        //shc_log_info(shc_logname, "read id %u, read len %u, start %u, seqlen %u\n",
        //                             ri_it->start, acc.len, ri_it->start, graph[vd].seq_len());


        char * node_left = acc.read_ptr + ri_it->start - 1;
        //shc_log_info(shc_logname, "L %c\n", *node_left);
        in_base_set.insert(*node_left);  //index start at 0
        char * node_right = acc.read_ptr + ri_it->start + curr_node.seq_len();
        out_base_set.insert(*node_right);  //index start at 0



        //shc_log_info(shc_logname, "R %c\n", *node_right);
        //shc_log_info(shc_logname, "***********one read bridge pass\n");
    }

    int bridged_in = boost::in_degree(vd, graph) - in_base_set.size();
    int bridged_out = boost::out_degree(vd, graph) - out_base_set.size();
    /*
    for(std::set<char>::iterator it=in_base_set.begin(); it!=in_base_set.end() ;it++)
    {
        shc_log_info(shc_logname, "in set %c\n", *it);;
    }
    for(std::set<char>::iterator it=out_base_set.begin(); it!=out_base_set.end() ;it++)
    {
        shc_log_info(shc_logname, "out set %c\n", *it);;
    }
    */

    in_base_set.clear();
    out_base_set.clear();

    if(bridged_in == 0 && bridged_out == 0)
    {
        //shc_log_info(shc_logname, "Finish bridged\n");
        return true;
    }
    else if ( bridged_in==1 && bridged_out==1 )
    {
        //shc_log_info(shc_logname, "Finish bridged\n");
        return true;
    }
    //shc_log_info(shc_logname, "Finish not bridged\n");
    return false;
}


void Sequence_graph_handler::bridge_all_xnodes()
{
    start_timer(&timer);
    shc_log_info(shc_logname, "Start bridging xnode\n");
    Node_set_iterator x_iter;
    int num_xnode_bridged = 0;
    for(x_iter=xnode_set.begin(); x_iter!=xnode_set.end(); ++x_iter)
    {
        refresh_bridging_reads(*x_iter);
        //log_node_info(*x_iter, false ,false);
        if(is_xnode_bridged(*x_iter))
        {
            //shc_log_info(shc_logname, "is bridging xnode\n");
            perform_bridging(*x_iter);
            num_xnode_bridged++;
        }
        else
        {
            //shc_log_info(shc_logname, "not bridging xnode\n");
        }
    }
    shc_log_info(shc_logname, "%d xnode bridged\n", num_xnode_bridged);
    shc_log_info(shc_logname, "Finish bridging xnode\n");
#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "multibridge finish ";
    stop_timer(&timer);
#endif
}

void Sequence_graph_handler::perform_bridging(vd_t vd)
{
    bundled_node_p & curr_node = graph[vd];
    assert(curr_node.reads_info.size() > 0);
    assert(boost::in_degree(vd, graph)>=2 && boost::out_degree(vd, graph)>=2);

    std::vector<vd_t> u_list, w_list;
    std::vector<vd_t>::iterator node_it;

    bool contain_loop = false;
    vd_t u_loop_node, w_loop_node;
    kmer_count_t loop_weight;

    // create left extension node
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "create left extension\n");
#endif
    in_eip_t vd_in_eip = boost::in_edges(vd, graph);
    for(in_ei_t it=vd_in_eip.first ;it!=vd_in_eip.second; it++)
    {
        vd_t u = create_left_node(*it);
        u_list.push_back(u);
        if(boost::source(*it, graph)==boost::target(*it, graph))
        {
            contain_loop = true;
            u_loop_node = u;
            loop_weight = graph[*it].weight;
        }
    }

    // create right extension node
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "create right extension\n");
#endif
    out_eip_t vd_out_eip = boost::out_edges(vd, graph);
    for(out_ei_t it=vd_out_eip.first ;it!=vd_out_eip.second; it++)
    {
        vd_t w = create_right_node(*it);
        w_list.push_back(w);
        if(boost::source(*it, graph)==boost::target(*it, graph))
            w_loop_node = w;
    }

#ifdef LOG_SEQ_GRAPH
    //log ulist wlist nodes
    shc_log_info(shc_logname, "ul%d\n", u_list.size());
    for(std::vector<vd_t>::iterator it=u_list.begin(); it!=u_list.end(); it++)
    {
        shc_log_info(shc_logname, "%s\n", graph[*it].seq.c_str());
    }

    shc_log_info(shc_logname, "wl%d\n", w_list.size());
    for(std::vector<vd_t>::iterator it=w_list.begin(); it!=w_list.end(); it++)
    {
        shc_log_info(shc_logname, "%s\n", graph[*it].seq.c_str());
    }
#endif


    // deal with self loop
    if (contain_loop)
    {
        shc_log_info(shc_logname, "contain self loop\n");
        aer_t aer = boost::add_edge(u_loop_node, w_loop_node, graph);
        graph[aer.first] = bundled_edge_p(loop_weight+2);
    }

    // clear edges around v, and vertex will be destroyed at the end
    boost::clear_vertex(vd, graph);

    // link extension nodes
    std::vector<Bdg_read_info> & reads_info = curr_node.reads_info;
    std::vector<Bdg_read_info>::iterator ri_it;

    std::vector<vd_t> matched_u;
    std::vector<vd_t> matched_w;

    std::set<vd_t> bridged_u;
    std::set<vd_t> bridged_w;

    // make sure reads of vertex are not repeated
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "bridging between extension nodes\n");
#endif
    for(ri_it=reads_info.begin(); ri_it!=reads_info.end(); ++ri_it)
    {
        Read_acc acc = read_list.get_read(ri_it->read_id);
        char * l_seq_start = acc.read_ptr + ri_it->start - 1;
        for(node_it=u_list.begin(); node_it!=u_list.end(); ++node_it)
        {
            bundled_node_p & node_p = graph[*node_it];
            //shc_log_info(shc_logname, "left ext %s\n", node_p.seq.c_str());
            if(memcmp(&node_p.seq.at(0), l_seq_start, node_p.seq_len())==0)
            {
               // shc_log_info(shc_logname, "matched left ext %s\n",node_p.seq.c_str());
                matched_u.push_back(*node_it);
                bridged_u.insert(*node_it);
            }
        }

        char * seq_start = acc.read_ptr + ri_it->start;
        for(node_it=w_list.begin(); node_it!=w_list.end(); ++node_it)
        {
            bundled_node_p & node_p = graph[*node_it];
            //shc_log_info(shc_logname, "right ext %s\n", node_p.seq.c_str());
            if(memcmp(&node_p.seq.at(0), seq_start, node_p.seq_len())==0)
            {
                //shc_log_info(shc_logname, "matched right ext %s\n", node_p.seq.c_str());
                matched_w.push_back(*node_it);
                bridged_w.insert(*node_it);
            }
        }
        if(matched_u.size()!=1 || matched_w.size()!=1)
        {
            shc_log_info(shc_logname, "in special case that some kmer in read is missing while buidling the graph\n");
            matched_u.clear();
            matched_w.clear();
            continue;
        }

        // create edge between nodes and link read
        vd_t u = matched_u[0];
        vd_t w = matched_w[0];
        node_link_read(u, ri_it->read_id, ri_it->start-1);
        node_link_read(w, ri_it->read_id, ri_it->start);

        //if there is no such edge
        aer_t aer = boost::edge(u,w,graph);
        if(!aer.second)
        {
            aer = boost::add_edge(u, w, graph);
            graph[aer.first] = bundled_edge_p(curr_node.seq_len());
        }
        matched_u.clear();
        matched_w.clear();
    }

    //log bridged nodes
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "check special case\n");
    shc_log_info(shc_logname, "ub%d\n", bridged_u.size());
    for(std::set<vd_t>::iterator it=bridged_u.begin(); it!=bridged_u.end(); it++)
    {
        shc_log_info(shc_logname, "%s\n", graph[*it].seq.c_str());
    }

    shc_log_info(shc_logname, "wb%d\n", bridged_w.size());
    for(std::set<vd_t>::iterator it=bridged_w.begin(); it!=bridged_w.end(); it++)
    {
        shc_log_info(shc_logname, "%s\n", graph[*it].seq.c_str());
    }
#endif

    //if this is the special x-node case
    if(u_list.size()-bridged_u.size()==1 &&
       w_list.size()-bridged_w.size()==1)
    {
        std::vector<vd_t>::iterator u_it, w_it;

        for(u_it=u_list.begin(); u_it!=u_list.end(); ++u_it)
        {
            if(bridged_u.find(*u_it)==bridged_u.end())
                break;
        }
        for(w_it=w_list.begin(); w_it!=w_list.end(); ++w_it)
        {
            if(bridged_w.find(*w_it)==bridged_w.end())
                break;
        }

        aer_t aer = boost::add_edge(*u_it, *w_it, graph);
        graph[aer.first] = bundled_edge_p(curr_node.seq_len());
#ifdef LOG_SEQ_GRAPH
        shc_log_info(shc_logname, "special bridging apply\n");
#endif
    }
    else
    {
        assert(u_list.size()-bridged_u.size()==0 &&
               w_list.size()-bridged_w.size()==0);
#ifdef LOG_SEQ_GRAPH
        shc_log_info(shc_logname, "everything directly bridged\n");
#endif
    }

    //condense graph
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "local condense bridged nodes\n");
#endif
    local_condense_bridged(u_list, w_list);

    boost::remove_vertex(vd, graph);
}

// This function takes one argument, condense left or right if it can,
// otherwise return false
bool Sequence_graph_handler::
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
            new_node.reads_info = source_node.reads_info;
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

            bundled_node_p & new_node = graph[vd_new];
            bundled_node_p & target_node = graph[vd_target];
            new_node.count = (new_node.count + target_node.count)/2;
            // info is only updated in the left side
            bundled_edge_p & edge_p = graph[*(out_eip.first)];
            new_node.seq = new_node.seq + target_node.seq.substr(edge_p.weight) ;
        }
    }

    assert(boost::out_degree(vd_source, graph) == 1 &&
           boost::in_degree(vd_target, graph) == 1);

    return true;
}

void Sequence_graph_handler::
condense_left(vd_t vd_target, vd_t vd_new, std::set<vd_t> & vd_remove_set)
{
    vd_t vd_source;
    while (condense(vd_target, true, vd_new))
    {
        //get the source node and clear
        in_eip_t in_eip = boost::in_edges(vd_target, graph);
        vd_source = boost::source(*(in_eip.first), graph);

        boost::clear_in_edges(vd_target, graph);
        vd_remove_set.insert(vd_source);

        //prepare for next iteration
        vd_target = vd_source;
    }
    vd_t vd_leftmost = vd_target;
    //link left nodes to new nodes
    in_eip_t in_eip = boost::in_edges(vd_leftmost, graph);
    for(in_ei_t it = in_eip.first; it!=in_eip.second; it++)
    {
        aer_t aer = boost::add_edge(boost::source(*it, graph), vd_new, graph);
        graph[aer.first] = bundled_edge_p(graph[*it].weight);
    }
    boost::clear_in_edges(vd_leftmost, graph);
}


void Sequence_graph_handler::
condense_right(vd_t vd_source, vd_t vd_new ,std::set<vd_t> & vd_remove_set)
{
    vd_t vd_target;
    while (condense(vd_source, false, vd_new))
    {
        //get the source node and clear
        out_eip_t out_eip = boost::out_edges(vd_source, graph);
        vd_target = boost::target(*(out_eip.first), graph);

        boost::clear_out_edges(vd_source, graph);
        vd_remove_set.insert(vd_target);

        //prepare for next iteration
        vd_source = vd_target;
    }

    vd_t vd_rightmost = vd_source; // the node which stands for the right end
    //link right node
    out_eip_t out_eip = boost::out_edges(vd_rightmost, graph);
    for(out_ei_t it = out_eip.first; it!=out_eip.second; it++)
    {
        aer_t aer = boost::add_edge(vd_new, boost::target(*it, graph), graph);
        graph[aer.first] = bundled_edge_p(graph[*it].weight);
    }
    boost::clear_out_edges(vd_rightmost, graph);
}

bool Sequence_graph_handler::is_condensable(vd_t vd, bool is_left)
{
    vd_t vd_source, vd_target;
    if(is_left)
    {
        if (boost::in_degree(vd, graph) != 1)
            return false;
        vd_target = vd;
        in_eip_t in_eip = boost::in_edges(vd, graph);
        vd_source = boost::source(*(in_eip.first), graph);
        return (boost::out_degree(vd_source, graph) == 1);
    }
    else
    {
        if (boost::out_degree(vd, graph) != 1)
            return false;
        vd_source = vd;
        out_eip_t out_eip = boost::out_edges(vd, graph);
        vd_target = boost::target(*(out_eip.first), graph);
        return (boost::in_degree(vd_target, graph) == 1);
    }
}

// always remember to clear the vd_remove_set afterward
void Sequence_graph_handler::
local_condense(vd_t vd, std::set<vd_t> & vd_remove_set)
{
    bool left_condensable = is_condensable(vd, true);
    bool right_condensable = is_condensable(vd, false);

    if( left_condensable || right_condensable)
    {
        vd_t vd_new = boost::add_vertex(graph);
        graph[vd_new] = graph[vd];
        graph[vd_new].reads_info = graph[vd].reads_info;
        //log_node_info(*n_it, true, false);

        std::reverse(graph[vd_new].seq.begin(), graph[vd_new].seq.end());
        condense_left(vd, vd_new, vd_remove_set);
        std::reverse(graph[vd_new].seq.begin(), graph[vd_new].seq.end());

        condense_right(vd, vd_new ,vd_remove_set);

        refresh_bridging_reads(vd_new);
        boost::clear_vertex(vd, graph);
        vd_remove_set.insert(vd);
    }
    else
    {
        ;//shc_log_info(shc_logname, "condense nothing\n");
    }
}

// all nodes have been bridged
void Sequence_graph_handler::
local_condense_bridged(std::vector<vd_t> & u_list, std::vector<vd_t> & w_list)
{
    //shc_log_info(shc_logname, "Start local condense\n");
    std::set<vd_t> vd_remove_set;
    std::vector<vd_t>::iterator n_it;

    u_list.insert(u_list.end(), w_list.begin(), w_list.end());

    for(n_it=u_list.begin(); n_it!=u_list.end(); ++n_it)
    {
        local_condense(*n_it, vd_remove_set);
    }

    //shc_log_info(shc_logname, "remove nodes associated with local condense\n");
    for(std::set<vd_t>::iterator it=vd_remove_set.begin(); it != vd_remove_set.end(); it++)
    {
        //shc_log_info(shc_logname, "in %d, out %d\n", boost::in_degree(*it, graph),
        assert(boost::in_degree(*it, graph)==0 &&
                boost::out_degree(*it, graph)==0);
        boost::clear_vertex(*it,graph);
        boost::remove_vertex(*it,graph);
    }
}

bool Sequence_graph_handler::
is_read_bridge_node(read_num_t info_read_id, read_length_t info_start, vd_t vd)
{
    Read_acc acc = read_list.get_read(info_read_id);
    if (info_start==0)
        return false;
    if (acc.len <= info_start + graph[vd].seq_len())
        return false;
    return memcmp(&graph[vd].seq.at(0),
            acc.read_ptr+info_start, graph[vd].seq_len())==0;
}

Sequence_graph_handler::vd_t Sequence_graph_handler::
create_left_node(ed_t ed)
{
    bundled_node_p & source_node = graph[boost::source(ed, graph)];
    bundled_node_p & curr_node = graph[boost::target(ed, graph)];

    vd_t u = boost::add_vertex(graph);
    aer_t new_aer = boost::add_edge(boost::source(ed, graph), u, graph);
    graph[new_aer.first] = bundled_edge_p(graph[ed].weight + 1);

    int left_index = source_node.seq.size() - 1 - graph[ed].weight;
    kmer_count_t count = (source_node.count + curr_node.count)/2;
    std::string seq = source_node.seq.at(left_index) + curr_node.seq;
    graph[u] = bundled_node_p(seq, count);

    std::vector<Bdg_read_info>::iterator it;
    for(it=source_node.reads_info.begin(); it!=source_node.reads_info.end(); ++it)
    {
        if(is_read_bridge_node(it->read_id, it->start-1, u))
        {
            node_link_read(u, it->read_id, it->start-1);
        }
    }
    return u;
}

Sequence_graph_handler::vd_t Sequence_graph_handler::
create_right_node(ed_t ed)
{
    bundled_node_p & curr_node = graph[boost::source(ed, graph)];
    bundled_node_p & target_node = graph[boost::target(ed, graph)];

    vd_t w = boost::add_vertex(graph);
    aer_t new_aer = boost::add_edge(w, boost::target(ed, graph), graph);
    graph[new_aer.first] = bundled_edge_p(graph[ed].weight + 1);

    kmer_count_t count = (curr_node.count + target_node.count)/2;
    std::string seq = curr_node.seq + target_node.seq.at(graph[ed].weight);
    graph[w] = bundled_node_p(seq, count);

    std::vector<Bdg_read_info>::iterator it;
    for(it=target_node.reads_info.begin(); it!=target_node.reads_info.end(); ++it)
    {
        if(is_read_bridge_node(it->read_id, it->start+1, w))
        {
            node_link_read(w, it->read_id, it->start+1);
        }
    }
    return w;
}

void Sequence_graph_handler::find_known_path(int max_hop)
{
    start_timer(&timer);
    shc_log_info(shc_logname, "Start find known path\n");
    typedef std::multimap<std::string, Node_index> Kmer_ni_MMap;
    typedef std::multimap<std::string, Node_index>::iterator
                                             Kmer_ni_MMap_iterator;
    typedef std::pair<Kmer_ni_MMap_iterator, Kmer_ni_MMap_iterator>
                                        Kmer_ni_MMap_iterator_pair;
    Kmer_ni_MMap kmer_ni_mmap;

    std::string kmer_base;
    kmer_base.resize(kmer_length);

    // construct kmer to node-index multimap, it is needed
    // since we perform multibridging, there can be multiple vertex that
    // that contain the start of a read if it is not a bridging read
    shc_log_info(shc_logname, "Building kmer->vd start_index multimap\n");
    vip_t vip = boost::vertices(graph);
    for(vi_t vi=vip.first; vi!=vip.second; ++vi)
    {
        bundled_node_p & curr_node = graph[*vi];
        for(int i=0; i<curr_node.seq_len()-kmer_length+1; ++i)
        {
            memcpy(&kmer_base.at(0), &curr_node.seq.at(i), kmer_length);
            kmer_ni_mmap.insert(std::make_pair(kmer_base, Node_index(*vi, i)));
        }
    }

    std::string k_start, k_end;
    k_start.resize(kmer_length);
    k_end.resize(kmer_length);
    int num_known_path = 0;
    //iterate through all read to get all known path
    //shc_log_info(shc_logname, "Iterating reads to find path\n");
    for(read_num_t i=0; i<read_list.get_num_read(); i++)
    {
        //shc_log_info(shc_logname, "start new read\n");
        Read_acc acc = read_list.get_read(i);
        memcpy(&k_start.at(0), acc.read_ptr, kmer_length);
        memcpy(&k_end.at(0),
                acc.read_ptr+acc.len-kmer_length, kmer_length);

        //if the front and end kmers in the read does not corresponds to nodes
        Kmer_ni_MMap_iterator_pair se_it = kmer_ni_mmap.equal_range(k_start);
        if(se_it.first==se_it.second)
            continue;
        se_it = kmer_ni_mmap.equal_range(k_end);
        if(se_it.first==se_it.second)
            continue;

        int curr_hop = max_hop;

        // find the corresponding start node
        Kmer_ni_MMap_iterator_pair kni_it_pair =
                                        kmer_ni_mmap.equal_range(k_start);
        //shc_log_info(shc_logname, "search matched nodes\n");
        for(Kmer_ni_MMap_iterator kni_it =kni_it_pair.first;
                                  kni_it!=kni_it_pair.second; ++kni_it)
        {
            vd_t vd = kni_it->second.vd;
            read_length_t val_node_index = kni_it->second.start;
            bundled_node_p & start_node = graph[vd];
            char * val_node_ptr = &start_node.seq.at(val_node_index);
            read_length_t val_node_len = start_node.seq_len() - val_node_index;

#ifdef LOG_SEQ_GRAPH
            std::string read_seq;
            read_seq.resize(acc.len);
            memcpy(&read_seq.at(0), acc.read_ptr, acc.len);
            shc_log_info(shc_logname, "Seq %s\n", read_seq.c_str());
#endif
            // if the first node matches the start of the read
            if(is_partial_match(acc.read_ptr, acc.len,
                                val_node_ptr, val_node_len))
            {
                std::vector<vd_t> path;
                if(search_seq(acc.read_ptr,
                              acc.len,
                              kni_it->second.vd,
                              kni_it->second.start, curr_hop, path))
                {
                    // there is a full cover path
                    //shc_log_info(shc_logname, "find matched path\n");
                    read_known_path_map[i] = path;
                    //known_path_set.insert(path);
                    num_known_path++;
                }
            }
        }
    }
    shc_log_info(shc_logname, "Finding %d full matching path\n", num_known_path);
    shc_log_info(shc_logname, "Finish find known path\n");
#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "finding " << num_known_path << " known path finish, ";
    stop_timer(&timer);
#endif
}


bool Sequence_graph_handler::
search_seq(const char * curr_seq_ptr, read_length_t curr_seq_len, vd_t vd,
                   read_length_t val_start, int curr_hop, std::vector<vd_t> & path)
{
#ifdef LOG_SEQ_GRAPH
    std::string path_seq;
    read_length_t seq_full_len = curr_seq_len;
#endif
    while((curr_hop--)>0)
    {
        path.push_back(vd);
        bundled_node_p & curr_node = graph[vd];
        assert(curr_node.seq_len() >= val_start);
        read_length_t val_node_len = curr_node.seq_len() - val_start;
#ifdef LOG_SEQ_GRAPH
        if(val_node_len!=0)
            path_seq.insert(path_seq.end(), curr_node.seq.begin()+val_start, curr_node.seq.end());
#endif
        if(curr_seq_len <= val_node_len)
        {
#ifdef LOG_SEQ_GRAPH
            shc_log_info(shc_logname, "success %u, %u\n", curr_seq_len, val_node_len);
            path_seq.resize(seq_full_len);
            shc_log_info(shc_logname, "         %s\n", path_seq.c_str());
#endif
            return true;
        }

        //shc_log_info(shc_logname, "seqlen %u, val node len %u\n", curr_seq_len, val_node_len);

        curr_seq_ptr += val_node_len;
        curr_seq_len -= val_node_len;

        // iterate out edges to find next match node
        bool find_match = false;
        out_eip_t out_eip = boost::out_edges(vd, graph);

        //std::string seq_debug;
        //seq_debug.resize(curr_seq_len);
        //memcpy(&seq_debug.at(0), curr_seq_ptr, curr_seq_len);
        //shc_log_info(shc_logname, "Seq %s\n", seq_debug.c_str());
        for(out_ei_t it = out_eip.first; it!=out_eip.second; it++)
        {
            read_length_t overlap = graph[*it].weight;
            vd_t vd_target = boost::target(*it, graph);
            bundled_node_p & next_node = graph[vd_target];

            //log_seq_only(vd_target, overlap);
            //the case when the current node seq contains the next seq
            if(next_node.seq_len() == overlap)
            {
                vd = vd_target;
                val_start = overlap;
                find_match = true;
                break;
            }

            char * val_node_ptr = &next_node.seq.at(overlap);

            read_length_t val_node_len = next_node.seq_len() - overlap;

            //shc_log_info(shc_logname, "val_node_len %d\n", val_node_len);

            if(is_partial_match(curr_seq_ptr, curr_seq_len,
                                val_node_ptr, val_node_len ))
            {
                vd = vd_target;
                val_start = overlap;
                find_match = true;
                break;
            }
        }
        //check if out nodes contain matched one
        if(!find_match)
        {
            //shc_log_info(shc_logname, "Read fails\n");
            return false;
        }
    }
    return false;
}

bool Sequence_graph_handler::
is_partial_match(const char * ptr1, const read_length_t len1,
              const char * ptr2, const read_length_t len2)
{
    return memcmp(ptr1, ptr2, std::min(len1, len2))==0;
}

void Sequence_graph_handler::break_self_loops()
{
    start_timer(&timer);
    shc_log_info(shc_logname, "break self loops\n");
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    vi_t next, it=vip.first;
    for(next=it; it!=vip.second; it=next)
    {
        next++;
        vd_t vd = *it;
        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            if(vd_target == vd)  // a self loop
            {
                shc_log_info("find a self loops, delete\n");
                boost::remove_edge(*out_ei, graph);
                break; //since for each node, it only contains one self node
            }
        }
    }
    shc_log_info(shc_logname, "Finish break self loops\n");
#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "break all self loop finish, ";
    stop_timer(&timer);
#endif
}

void Sequence_graph_handler::break_all_cycles()
{
    start_timer(&timer);
    shc_log_info(shc_logname, "Start break all cycle\n");

    shc_log_info(shc_logname,"current number of vertices is %d\n",
                                                boost::num_vertices(graph));
    std::set<vd_t> acyclic_node_set;
    std::deque<vd_t> cycle_path;
    int i=0;
    while (find_cycle(acyclic_node_set, cycle_path))
    {
        //shc_log_info(shc_logname, "Start break cycle with size %d\n", cycle_path.size());
        break_cycle(cycle_path);
        cycle_path.clear();
        //std::cout << "finish break " << (++i) << " cycles" << std::endl;
    }
    condense_graph();

    shc_log_info(shc_logname, "Finish break all cycle\n");
#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "break all cycles finish, ";
    stop_timer(&timer);
#endif
}

bool Sequence_graph_handler::
find_cycle(std::set<vd_t> & acyclic_node_set, std::deque<vd_t> & cycle_path)
{
    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        std::set<vd_t>::iterator an_it = acyclic_node_set.find(vd);
        if(an_it == acyclic_node_set.end()) //not in the set
        {
            if(is_node_inside_cycle(vd, acyclic_node_set, cycle_path))
            {
                //if there is a cycle
                 return true;
            }
        }
    }
    return false;
}

bool Sequence_graph_handler::
is_node_inside_cycle(vd_t vd, std::set<vd_t> & acyclic_node_set, std::deque<vd_t> & cycle_path)
{
    //shc_log_info(shc_logname, "\t\tchecking a new node\n");
    cycle_path.push_back(vd);
    //shc_log_info(shc_logname, "pushed a node %s\n", graph[vd].seq.c_str());
    // for all out edges
    out_eip_t out_eip = boost::out_edges(vd, graph);
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        std::deque<vd_t>::iterator cycle_index_it;
        cycle_index_it = std::find(cycle_path.begin(), cycle_path.end(), vd_target);
        //if find same vertex in the middle
        if(cycle_index_it  != cycle_path.end())
        {
            cycle_path.push_back(vd_target);
            //shc_log_info(shc_logname, "cycle node %s\n", graph[vd_target].seq.c_str());
            if(*cycle_index_it != *(cycle_path.begin()) )
            {
                cycle_path.erase(cycle_path.begin(), cycle_index_it);
            }
            //assert((cycle_path.back()) == (cycle_path.front()));
            return true;
        }
        if(acyclic_node_set.find(vd_target) != acyclic_node_set.end())
            continue;

        return is_node_inside_cycle(vd_target, acyclic_node_set, cycle_path);
    }
    acyclic_node_set.insert(vd);
    cycle_path.clear();
    //shc_log_info(shc_logname, "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n");
    return false;
}

void Sequence_graph_handler::break_cycle(std::deque<vd_t> & cycle_path)
{
    // find anything to break
    cycle_path.pop_front();
    std::vector<vd_t> xnodes;
    std::vector<vd_t> prev_node;
    std::vector<vd_t> next_node;
    for(std::deque<vd_t>::iterator c_it=cycle_path.begin();
                                    c_it!=cycle_path.end(); ++c_it)
    {
        vd_t vd = *c_it;
        bundled_node_p & curr_node = graph[vd];
        if(boost::in_degree(vd, graph)>=2 && boost::out_degree(vd, graph)>=2)
        {
            xnodes.push_back(vd);
            if(c_it==cycle_path.begin())
                prev_node.push_back(*(cycle_path.end()-1));
            else
                prev_node.push_back(*(c_it-1));
            if(c_it==cycle_path.end()-1)
                next_node.push_back(*(cycle_path.begin()));
            else
                next_node.push_back(*(c_it+1));
        }
    }
    if(xnodes.size()==0)
    {
        for(std::deque<vd_t>::iterator c_it=cycle_path.begin();
                                        c_it!=cycle_path.end(); ++c_it)
        {
            vd_t vd = *c_it;
            bundled_node_p & curr_node = graph[vd];
            if(boost::in_degree(vd, graph)>=2 || boost::out_degree(vd, graph)>=2)
            {
                xnodes.push_back(vd);
                if(c_it==cycle_path.begin())
                    prev_node.push_back(*(cycle_path.end()-1));
                else
                    prev_node.push_back(*(c_it-1));
                if(c_it==cycle_path.end()-1)
                    next_node.push_back(*(cycle_path.begin()));
                else
                    next_node.push_back(*(c_it+1));
            }
        }
    }
    if(xnodes.size()==0)
    {
        shc_log_warning("no enough xnodes in cycle, unable to break\n");
        exit(1);
    }
    // finding where to break
    int min_disrupt = std::numeric_limits<int>::max();
    int sum_copy_count = 0;
    int sum_c_in_and_o_out = 0;
    vd_t vd_split, vd_c_in, vd_c_out;   // those are nodes inside cycle

    //shc_log_info(shc_logname, "xnodes size %d, cycle size %d\n", xnodes.size(), cycle_path.size());
    for(std::vector<vd_t>::iterator c_it=xnodes.begin(); c_it!=xnodes.end(); ++c_it)
    {
        vd_t vd_curr = *c_it;
        vd_t vd_prev = prev_node.at(c_it-xnodes.begin());
        vd_t vd_next = next_node.at(c_it-xnodes.begin());

        int o_in_copy_count_sum=0, c_in_copy_count=0;
        in_eip_t in_eip = boost::in_edges(vd_curr, graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
            vd_t vd_source = boost::source(*in_ei, graph);
            if(vd_source != vd_prev)
                o_in_copy_count_sum += graph[*in_ei].copy_count;
            else
                c_in_copy_count = graph[*in_ei].copy_count;
        }

        int  o_out_copy_count_sum=0, c_out_copy_count=0;
        out_eip_t out_eip = boost::out_edges(vd_curr, graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            if(vd_target != vd_next)
                o_out_copy_count_sum += graph[*out_ei].copy_count;
            else
                c_out_copy_count = graph[*out_ei].copy_count;
        }

        int disrupt = std::abs(o_in_copy_count_sum-c_out_copy_count) +
                      std::abs(o_out_copy_count_sum-c_in_copy_count);
        //shc_log_info(shc_logname, "disrupt %d\n", disrupt);
        if(min_disrupt > disrupt)
        {
            min_disrupt = disrupt;
            vd_split = vd_curr;
            vd_c_in = vd_prev;
            vd_c_out = vd_next;
            sum_copy_count = o_in_copy_count_sum + o_out_copy_count_sum +
                             c_out_copy_count + c_in_copy_count;
            sum_c_in_and_o_out = c_in_copy_count + o_out_copy_count_sum;
        }
    }

    break_cycle_splitting(vd_split, vd_c_in, vd_c_out, sum_copy_count,
                                                        sum_c_in_and_o_out);
}

void Sequence_graph_handler::
break_cycle_splitting(vd_t vd_split, vd_t vd_c_in, vd_t vd_c_out,
                                    int sum_copy_count, int sum_c_in_and_o_out)
{
    // actual breaking
    //shc_log_info(shc_logname, "breaking cycle node %s\n", graph[vd_split].seq.c_str());
    vd_t vd_clone = boost::add_vertex(graph);
    double proportion = 0.5;
    if ( sum_copy_count != 0 )
    {
        proportion = static_cast<double>(sum_c_in_and_o_out) /
                     static_cast<double>(sum_copy_count);
    }
    kmer_count_t clone_count =
                    static_cast<kmer_count_t>(proportion*graph[vd_split].count);

    graph[vd_clone] = bundled_node_p(graph[vd_split].seq, clone_count);
    graph[vd_clone].reads_info = graph[vd_split].reads_info;
    graph[vd_split].count =
                static_cast<kmer_count_t>((1-proportion)*graph[vd_split].count);

    aer_t aer_cycle_in = boost::edge(vd_c_in, vd_split, graph);
    assert(aer_cycle_in.second);
    aer_t aer_cycle_out = boost::edge(vd_split, vd_c_out, graph);
    assert(aer_cycle_out.second);

    aer_t aer_new = boost::add_edge(vd_c_in, vd_clone, graph);
    graph[aer_new.first] = bundled_edge_p(graph[aer_cycle_in.first].weight);
    boost::remove_edge(aer_cycle_in.first, graph);

    out_eip_t out_eip = boost::out_edges(vd_split, graph);
    std::vector<ed_t> edge_to_remove;
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        if(vd_target != vd_c_out)
        {
            aer_new = boost::add_edge(vd_clone, vd_target, graph);
            graph[aer_new.first] = bundled_edge_p(graph[aer_cycle_out.first].weight);
            edge_to_remove.push_back(*out_ei);
        }
    }

    for(std::vector<ed_t>::iterator it=edge_to_remove.begin();
                                    it!=edge_to_remove.end(); ++it)
    {
        boost::remove_edge(*it, graph);
    }
    //refresh_bridging_reads(vd_split);
    //refresh_bridging_reads(vd_clone);
    //shc_log_info(shc_logname, "split node finish\n");
}

void Sequence_graph_handler::find_mate_pairs()
{
    for(Read_path_map_iterator it = read_known_path_map.begin();
                               it !=read_known_path_map.end(); ++it)
    {
        ;
    }
}


void Sequence_graph_handler::log_all_node(bool is_log_node, bool is_log_read)
{
    shc_log_info(shc_logname, "log all nodes\n");
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        log_node_info(*it, is_log_node, is_log_read);
    }
}

void Sequence_graph_handler::log_xnode_set(bool is_log_node, bool is_log_read)
{
    shc_log_info(shc_logname, "log all xnodes\n");
    for(std::set<vd_t>::iterator it=xnode_set.begin(); it!=xnode_set.end(); it++)
    {
        log_node_info(*it, is_log_node, is_log_read);
    }
}

void Sequence_graph_handler::log_read_seq(Read_acc acc)
{
    std::string read;
    read.resize(acc.len);
    memcpy(&read.at(0), acc.read_ptr, acc.len);
    shc_log_info(shc_logname, "log read: %s\n", read.c_str());
}

void Sequence_graph_handler::log_seq_only(vd_t vd, read_length_t start)
{
    bundled_node_p & curr_node = graph[vd];
    shc_log_info(shc_logname, "********************************************************\n");
    shc_log_info(shc_logname, "Node start %d: %s\n", start, curr_node.seq.c_str());
}

void Sequence_graph_handler::log_node_info(vd_t vd, bool is_log_nodes ,bool is_log_read_info)
{
    bundled_node_p & curr_node = graph[vd];

    shc_log_info(shc_logname, "********************************************************\n");
    shc_log_info(shc_logname, "Node info: %s\n", curr_node.seq.c_str());
    shc_log_info(shc_logname, "in degree %d, out degree %d\n",
            boost::in_degree(vd,graph), boost::out_degree(vd,graph));
    if(is_log_nodes)
    {
        shc_log_info(shc_logname, "In nodes:\n");
        in_eip_t in_eip = boost::in_edges(vd,graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
            vd_t vd_source = boost::source(*in_ei, graph);
            shc_log_info(shc_logname, "%d %d %d %s\n", boost::in_degree(vd_source, graph),
                    boost::out_degree(vd_source, graph), graph[*in_ei].weight, graph[vd_source].seq.c_str());
        }
        shc_log_info(shc_logname, "Out nodes:\n");
        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            shc_log_info(shc_logname, "%d %d %d %s\n", boost::in_degree(vd_target, graph),
                    boost::out_degree(vd_target, graph) , graph[*out_ei].weight ,graph[vd_target].seq.c_str());
        }
    }

    shc_log_info(shc_logname, "%d Read cover it\n", curr_node.reads_info.size());
    if(is_log_read_info)
    {
        std::vector<Bdg_read_info> & reads_info = curr_node.reads_info;
        for(std::vector<Bdg_read_info>::iterator it=reads_info.begin();
                it!=reads_info.end(); it++)
        {
            //shc_log_info(shc_logname, "total_num read %u\n", read_list.get_num_read());

            shc_log_info(shc_logname, "read id %u, start %u\n", it->read_id, it->start);
            log_read_seq(read_list.get_read(it->read_id));
        }
    }
    assert(boost::in_degree(vd,graph)<=4 && boost::out_degree(vd,graph)<=4);
}

void Sequence_graph_handler::log_edge_weight()
{
    eip_t eip = boost::edges(graph);
    for(ei_t ei=eip.first; ei!=eip.second; ei++)
    {
        shc_log_info(shc_logname, "weight %d\n", graph[*ei].weight);
    }
}
