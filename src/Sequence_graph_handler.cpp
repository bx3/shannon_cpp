#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <bits/stl_vector.h>
#include <unistd.h>

#include "Sequence_graph_handler.h"

void Sequence_graph_handler::run_it(int comp_i, bool is_single_component)
{
    is_run_single_component = is_single_component;
    is_multi_thread = true;

    //std::cout <<" \t\tworking on comp " << comp_i << std::endl;
    curr_comp  = comp_i;
    //curr_comp = comp_i;
    if (curr_comp >= 0)
    {
        setup_input_file_path(comp_i);
    }

    coll_read_list.set_comp_id(comp_i); // for logging

    Block_timer main_timer;
    Block_timer dump_timer;
    start_timer(&main_timer);
        //build_kmer_graph_from_reads();

    uint64_t build_time, build_num_node, build_num_edge;
    uint64_t condense_time, condense_num_node, condense_num_edge;
    uint64_t rm_sus_time, rm_sus_num_node, rm_sus_num_edge;
    uint64_t collapse_time, collapse_num_node, collapse_num_edge;
    uint64_t update_xnode_set_time;
    uint64_t assign_read_time;
    uint64_t mb_time, mb_num_node, mb_num_edge;
    uint64_t cycle_break_time, cycle_break_num_node, cycle_break_num_edge;
    uint64_t find_known_path_time;
    uint64_t output_time;

    Block_timer part_timer;
    start_timer(&part_timer);

    build_kmer_graph_from_edges();
    build_num_node = get_num_nodes();
    build_num_edge = get_num_edges();
#ifdef PRINT_TIME
    std::cout << "build graph, load "<< coll_read_list.get_num_reads() << " read take, " << std::endl;
    stop_timer(&part_timer);
#endif
    stop_timer_np(&part_timer);
    build_time = part_timer.nTime;

    start_timer(&part_timer);
    condense_graph();
    condense_num_node = get_num_nodes();
    condense_num_edge = get_num_edges();
#ifdef PRINT_TIME
    std::cout << "condense graph take, " << std::endl;
    stop_timer(&part_timer);
#endif
    stop_timer_np(&part_timer);
    condense_time = part_timer.nTime;

#ifdef PRINT_GRAPH_STATS
    std::cout << "after condense num node " << (get_num_nodes()) << " num edge "
              << get_num_edges() << std::endl;
#endif

#ifdef PRINT_TIME
    start_timer(&part_timer);
#endif
    remove_all_suspicious_nodes();
    rm_sus_num_node = get_num_nodes();
    rm_sus_num_edge = get_num_edges();
#ifdef PRINT_TIME
    std::cout << "remove_all_suspicious_nodes, " << std::endl;
    stop_timer(&part_timer);
#endif
    stop_timer_np(&part_timer);
    rm_sus_time = part_timer.nTime;


#ifdef PRINT_GRAPH_STATS
    std::cout << "after rem suspicuous node " << (get_num_nodes()) << " num edge "
                      << get_num_edges() << std::endl;
#endif

// start
    start_timer(&part_timer);
    collapse_all();
    collapse_num_node = get_num_nodes();
    collapse_num_edge = get_num_edges();
#ifdef PRINT_TIME
    std::cout << "collapse_all, " << std::endl;
    stop_timer(&part_timer);
#endif
    stop_timer_np(&part_timer);
    collapse_time = part_timer.nTime;

#ifdef PRINT_GRAPH_STATS
    std::cout << "after collapse node " << (get_num_nodes()) << " num edge "
                              << get_num_edges() << std::endl;
#endif

    start_timer(&part_timer);
    update_xnode_set();
    int num_xnode = xnode_set.size();
#ifdef PRINT_TIME
    std::cout << "update_xnode_set, " << std::endl;
    stop_timer(&part_timer);
#endif
    stop_timer_np(&part_timer);
    update_xnode_set_time = part_timer.nTime;


    start_timer(&part_timer);
    assign_read_to_xnode();

#ifdef PRINT_TIME
    std::cout << "assign_read_to_xnode, " << std::endl;
    stop_timer(&part_timer);
#endif
    stop_timer_np(&part_timer);
    assign_read_time = part_timer.nTime;
    //for(std::set<vd_t>::iterator it=xnode_set.begin();
    //                    it!=xnode_set.end(); it++)
    //{
    //    std::cout << "node has num read " << graph[*it].reads_info.size() << std::endl;
    //}


    start_timer(&part_timer);
    bridge_all_xnodes();
    int num_xnode_remain = get_num_xnodes();
    //std::cout << "num_xnode_remain " << num_xnode_remain << std::endl;
    mb_num_node = get_num_nodes();
    mb_num_edge = get_num_edges();
#ifdef PRINT_TIME
    std::cout << "bridge_all_xnodes, " << std::endl;
    stop_timer(&part_timer);
#endif
    stop_timer_np(&part_timer);
    mb_time = part_timer.nTime;

#ifdef PRINT_GRAPH_STATS
    std::cout << "after bridge num node " << (get_num_nodes()) << " num edge "
              << get_num_edges() << std::endl;
#endif

    condense_graph();

    find_approximate_copy_count();

    start_timer(&part_timer);
    //break_self_loops();

#ifdef PRINT_TIME
    std::cout << "break_self_loops, " << (get_num_nodes()) << " num edge "
              << get_num_edges() << std::endl;
    stop_timer(&part_timer);
#endif
    //condense_graph();
    //std::cout << "after condense " << (get_num_nodes()) << " num edge "
    //          << get_num_edges() << std::endl;

    //log_graph_to_file(graph_writer , 102);
/*
    special_condense_self_loop();
    condense_graph();

#ifdef PRINT_GRAPH_STATS
    std::cout << "after special condense node " << (get_num_nodes()) << " num edge "
              << get_num_edges() << std::endl;
#endif
*/
    //log_all_node(true, false);



    start_timer(&part_timer);
    //py_break_cycles();
    break_all_cycles();

    cycle_break_num_node = get_num_nodes();
    cycle_break_num_edge = get_num_edges();
#ifdef PRINT_TIME
    std::cout << "break_all_cycles, " << (get_num_nodes()) << " num edge "
              << get_num_edges() << std::endl;
    stop_timer(&part_timer);
#endif
    stop_timer_np(&part_timer);
    cycle_break_time = part_timer.nTime;

    //int before_rm_node_num = get_num_nodes();
    //int before_rm_edge_num = get_num_edges();

    //remove_nodes_edges_if_not_cover_by_reads();

    //int after_rm_node_num = get_num_nodes();
    //int after_rm_edge_num = get_num_edges();

    //std::cout << "before_rm_node_num " << before_rm_node_num << std::endl;
    //std::cout << "before_rm_edge_num " << before_rm_edge_num << std::endl;
    //std::cout << "after_rm_node_num " << after_rm_node_num << std::endl;
    //std::cout << "after_rm_edge_num " << after_rm_edge_num << std::endl;

    find_approximate_copy_count();

    start_timer(&part_timer);
    //log_graph_to_file(graph_writer , 103);


    int num_known_path_found = find_known_path(max_hop_path);
#ifdef PRINT_TIME
    std::cout << "find_known_path, " << std::endl;
    stop_timer(&part_timer);
#endif
    stop_timer_np(&part_timer);
    find_known_path_time = part_timer.nTime;

    find_copy_count();

    Local_files & lf = setting.local_files;

    std::string node_dir =
        lf.output_seq_graph_path + lf.node_prefix + std::to_string(comp_i);
    std::string edge_dir =
        lf.output_seq_graph_path + lf.edge_prefix + std::to_string(comp_i);
    std::string path_dir =
        lf.output_seq_graph_path + lf.path_prefix + std::to_string(comp_i);

    //std::cout << "output node file " << node_dir << std::endl;
    //std::cout << "output edge file " << edge_dir << std::endl;
    //std::cout << "output path file " << path_dir << std::endl;

    //std::cout << "num isolated nodes " << get_num_isolated_nodes() << std::endl;

    start_timer(&dump_timer);
    //log_classify_edge_types();
    //log_graph_to_file(graph_writer , 104);
    output_components(node_dir, edge_dir, path_dir);
#ifdef PRINT_TIME
    std::cout << "output_components, " << std::endl;
    stop_timer(&dump_timer);
#endif
    stop_timer_np(&dump_timer);
    output_time = dump_timer.nTime;
    //clear_seq_graph_mem();

    if(is_single_component)
    {
        std::cout << std::endl;
        std::cout << "Finish component " << comp_i << ", SUMMARY " << std::endl;
        std::cout << "build_time            " << build_time<< " sec, Node: "
                  << build_num_node << " Edge: " << build_num_edge  << std::endl;
        std::cout << "condense_time         " << condense_time << " sec, Node: "
                  << condense_num_node << " Edge: " << condense_num_edge  << std::endl;
        std::cout << "rm_sus_time           " << rm_sus_time << " sec, Node: "
                  << rm_sus_num_node << " Edge: " << rm_sus_num_edge  << std::endl;
        std::cout << "collapse_time         " << collapse_time << " sec, Node: "
                  << collapse_num_node << " Edge: " << collapse_num_edge  << std::endl;
        std::cout << "update_xnode_set_time " << update_xnode_set_time << " sec, num xnode " << num_xnode << std::endl;
        std::cout << "assign_read_time      " << assign_read_time << " sec, No. assigned reads " << num_assigned_read << " out of total reads " << coll_read_list.get_num_reads() << std::endl;
        std::cout << "resolve_pair_time     " << resolve_pair_timer.nTime << " sec" << std::endl;
        std::cout << "mb_time               " << mb_time << " sec, Node: "
                  << mb_num_node << " Edge: " << mb_num_edge << ", "
                  << total_num_bridged << " nodes bridged"   << std::endl;
        std::cout << "cycle_break_time      " << cycle_break_time << " sec, Node: "
                  << cycle_break_num_node << " Edge: " << cycle_break_num_edge  << std::endl;
        std::cout << "find_known_path_time  " << find_known_path_time << " sec, num found " << num_known_path_found << std::endl;
        std::cout << "output_time           " << output_time << " sec, "
                  << num_single_seq << " single seq" <<  std::endl;
        std::cout << "Total time, ";
        stop_timer(&main_timer);
    }
    shc_log_info(shc_logname, "Finish comp %d, SUMMARY\n", comp_i);
    shc_log_info(shc_logname, "build_time            %d sec, Node: %u Edges: %u\n",
                    build_time, build_num_node, build_num_edge);
    shc_log_info(shc_logname, "condense_time         %d sec, Node: %u Edges: %u\n",
                    condense_time, condense_num_node, condense_num_edge);
    shc_log_info(shc_logname, "rm_sus_time           %d sec, Node: %u Edges: %u\n",
                    rm_sus_time, rm_sus_num_node, rm_sus_num_edge);
    shc_log_info(shc_logname, "collapse_time         %d sec, Node: %u Edges: %u\n",
                    collapse_time, collapse_num_node, collapse_num_edge);
    shc_log_info(shc_logname, "update_xnode_set_time %d sec\n", update_xnode_set_time);
    shc_log_info(shc_logname, "assign_read_time      %d sec, %d reads\n", assign_read_time, coll_read_list.get_num_reads());
    shc_log_info(shc_logname, "resolve pair time %d sec\n", resolve_pair_timer.nTime);
    shc_log_info(shc_logname, "mb_time               %d sec, Node: %u Edges: %u, %u node bridged, \n",
                    mb_time, mb_num_node, mb_num_edge, total_num_bridged);
    shc_log_info(shc_logname, "cycle_break_time      %d sec, Node: %u Edges: %u\n",
                    cycle_break_time, cycle_break_num_node, cycle_break_num_edge);
    shc_log_info(shc_logname, "find_known_path_time  %d sec num found %d\n", find_known_path_time, num_known_path_found);
    shc_log_info(shc_logname, "output_time           %d sec, %u single seq\n", output_time, num_single_seq);
    shc_log_info(shc_logname, "total time \n");
    log_stop_timer(&main_timer);
}

void Sequence_graph_handler::
setup_input_file_path(std::string kmer_path_, std::string s_read_path_,
                        std::string p1_read_path_, std::string p2_read_path_)
{
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start Custom file setup\n");
#endif
    glob_node_id = NODE_ID_NORMAL_START;
    kmer_path = kmer_path_;

    if(!s_read_path_.empty())
    {
        read_path_single_prefix = s_read_path_;
        coll_read_list.s_reads.setup(
            get_num_seq(s_read_path_));
    }

    if(!p1_read_path_.empty() && !p2_read_path_.empty())
    {
        read_path_p1_prefix = p1_read_path_;
        coll_read_list.p1_reads.setup(
            get_num_seq(p1_read_path_));
        read_path_p2_prefix = p2_read_path_;
        coll_read_list.p2_reads.setup(
            get_num_seq(p2_read_path_));
    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Finish Custom file setup\n");
#endif
}

size_t Sequence_graph_handler::
get_total_num_reads(std::string file_path_prefix)
{
    int i=0;
    size_t total_num_line = 0;
    std::string file_path = file_path_prefix + "_" + std::to_string(i);

    while(exist_path(file_path))
    {
        total_num_line += get_num_seq(file_path);
        i++;
        file_path = file_path_prefix + "_" + std::to_string(i);
    }
    //assert(total_num_line!=0);
    return total_num_line;
}

void Sequence_graph_handler::setup_input_file_path(int comp_i)
{
    kmer_path = setting.local_files.output_components_kmer_dir
              + setting.local_files.comp_kmer_prefix
              + std::to_string(comp_i);
    //std::cout << "kmer_path " << kmer_path << std::endl;
    if(has_single)
    {
        read_path_single_prefix = setting.local_files.output_components_read_dir
                         + setting.local_files.comp_read_prefix
                         + std::to_string(comp_i);
        size_t file_size =
            get_total_num_reads(read_path_single_prefix);
        //std::cout << "comp " << comp_i << "single file size " << file_size << std::endl;
        coll_read_list.s_reads.setup(file_size);
    }
    if(has_pair)
    {
        read_path_p1_prefix = setting.local_files.output_components_read_dir
                     + setting.local_files.comp_read_prefix
                     + std::to_string(comp_i)+"_p1";
        coll_read_list.p1_reads.setup(
            get_total_num_reads(read_path_p1_prefix));

        read_path_p2_prefix = setting.local_files.output_components_read_dir
                     + setting.local_files.comp_read_prefix
                     + std::to_string(comp_i)+"_p2";
        coll_read_list.p2_reads.setup(
            get_total_num_reads(read_path_p2_prefix));
    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "input path setup complete\n");
#endif
}

void Sequence_graph_handler::clear_seq_graph_mem()
{
    graph.clear();      // indeed free mem
    xnode_set.clear();
    paired_node_set.clear();
    p1_end.clear();
    //p1_end.shrink_to_fit();
    p2_front.clear();
    //p2_front.shrink_to_fit();
    known_path_set.clear();
    coll_read_list.clear();
}


void Sequence_graph_handler::build_kmer_graph_from_edges()
{
    std::ifstream file_reader(kmer_path.c_str());
    std::string count_str, kmer_base;
    double count;

    Kmer_Node_map kmer_node_map;
    kmer_node_map.set_empty_key("");

    // create nodes and
    while(  std::getline(file_reader, kmer_base, '\t') &&
            std::getline(file_reader, count_str)   )
    {
        count = std::stoi(count_str);
        kmer_length = kmer_base.size()-1;
        //std::cout << "kmer_length " << (uiunt16_t)kmer_length << std::endl;
        std::string base1 = kmer_base.substr(0, kmer_length);
        std::string base2 = kmer_base.substr(1, kmer_length);

        Kmer_Node_map_iterator kn_it_1 = kmer_node_map.find(base1);

        vd_t vd1,vd2;
        if (kn_it_1 == kmer_node_map.end())
        {
            vd1 = boost::add_vertex(graph);
            graph[vd1] = bundled_node_p(base1 ,glob_node_id++);
            kmer_node_map[base1] = vd1;
            exist_vd.insert(vd1);
        }
        else
        {
            vd1 = kn_it_1->second;
        }

        Kmer_Node_map_iterator kn_it_2 = kmer_node_map.find(base2);
        if (kn_it_2 == kmer_node_map.end())
        {
            vd2 = boost::add_vertex(graph);
            graph[vd2] = bundled_node_p(base2 ,glob_node_id++);
            kmer_node_map[base2] = vd2;
            exist_vd.insert(vd2);
        }
        else
        {
            vd2 = kn_it_2->second;
        }

        aer_t aer = boost::edge(vd1, vd2, graph);
        if (aer.second)
            graph[aer.first].count += count;
        else
        {
            aer = boost::add_edge(vd1, vd2, graph);
            graph[aer.first] = bundled_edge_p(kmer_length-1, count);
        }
    }

    // add node prevalence

    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        uint64_t out_sum_p = 0;
        out_eip_t out_eip = boost::out_edges(*it, graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            out_sum_p += graph[*out_ei].count;
        }

        uint64_t in_sum_p = 0;
        in_eip_t in_eip = boost::in_edges(*it ,graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
            in_sum_p += graph[*in_ei].count;
        }


        graph[*it].prevalence = (in_sum_p+out_sum_p)/2;
    }

    load_all_read(kmer_node_map);
}


vd_t Sequence_graph_handler::
process_one_read_to_build_graph( Kmer_Node_map& kmer_node_map,
                std::string & read_base, uint8_t node_type)
{
    std::string kmer_base;
    kmer_base.resize(kmer_length);
    vd_t prev_vd = NULL, vd_end = NULL, vd_front = NULL;
    bool is_front_node_set = false;
    bool is_kmer_missed = false;
    for(int i=0; i<read_base.size()-kmer_length+1; i++)
    {
        memcpy(&kmer_base.at(0), &read_base.at(i), kmer_length);
        Kmer_Node_map_iterator kn_iter = kmer_node_map.find(kmer_base);
        // If there is such node
        if(kn_iter != kmer_node_map.end())
        {
            vd_t & curr_vd = kn_iter->second;
            if(node_type == PAIR_END_NODE)
                vd_end = curr_vd;
            if(!is_front_node_set && node_type == PAIR_FRONT_NODE)
            {
                vd_front = curr_vd;
                is_front_node_set = true;
            }
            if (prev_vd != NULL)
            {
                aer_t aer = boost::edge(prev_vd, curr_vd, graph);
                if(!aer.second)
                {
                    aer = boost::add_edge(prev_vd, curr_vd, graph);
                    graph[aer.first] = bundled_edge_p(kmer_length-1, 1);
                }
                else
                {
                    // increment edge count
                    graph[aer.first].count++;
                }
            }
            prev_vd = curr_vd;
        }
        else
        {
            is_kmer_missed = true;
            prev_vd = NULL; //shc_log_info(shc_logname, "broken read\n");
        }
    }

    //if(read_base.size() > setting.single_read_length)
    //    shc_log_info(shc_logname, "%d: %s\n", read_base.size(), read_base.c_str());

    if(node_type == PAIR_END_NODE && !is_kmer_missed)
        return vd_end;
    else if(node_type == PAIR_FRONT_NODE && !is_kmer_missed)
        return vd_front;
    else
        return NULL;
}

void Sequence_graph_handler::load_all_read(Kmer_Node_map & kmer_node_map)
{
    if(has_single)
    {
        //std::cout << "read_path_single_prefix " << read_path_single_prefix << std::endl;
        if(exist_path(read_path_single_prefix)) //load direct path first
        {
            //std::cout << "load_all_single_read " << std::endl;
            load_all_single_read(read_path_single_prefix, kmer_node_map);
        }
        else
        {
            int i=0;
            std::string file_path = read_path_single_prefix + "_" + std::to_string(i);
            while(exist_path(file_path))
            {
                //std::cout << "load read " << read_path_single_prefix << std::endl;
                load_all_single_read(file_path, kmer_node_map);
                i++;
                file_path = read_path_single_prefix + "_" + std::to_string(i);

            }
        }
    }
    if(has_pair)
    {
        if(exist_path(read_path_p1_prefix) && exist_path(read_path_p2_prefix))
        {
            //std::cout << "load_all_paired_read_no_concat " << std::endl;
            load_all_paired_read_no_concat(read_path_p1_prefix,
                            read_path_p2_prefix, kmer_node_map);
        }
        else
        {
            load_all_paired_read(kmer_node_map);
        }
    }

    coll_read_list.declare_read_finish();
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "s_reads %u unique reads\n", coll_read_list.s_reads.get_num_reads());
    shc_log_info(shc_logname, "p1_reads %u unique reads\n", coll_read_list.p1_reads.get_num_reads());
    shc_log_info(shc_logname, "p2_reads %u unique reads\n", coll_read_list.p2_reads.get_num_reads());
    shc_log_info(shc_logname, "%u unique reads in read pool\n", coll_read_list.get_num_reads());
#endif
}

//return is valid
bool Sequence_graph_handler::trim_read(std::string & read, uint64_t i_read)
{
    // in case the next line char is readed
    if(read.back() != 'A' && read.back() != 'C' &&
       read.back() != 'G' && read.back() != 'T'  )
    {
        //shc_log_info(shc_logname, "B read %s\n", read.c_str());
        read.resize(read.size()-1);
        //shc_log_info(shc_logname, "A read %s\n", read.c_str());
    }
    boost::trim(read);
    for(int i=0; i<read.size(); i++)
    {
        char letter = read.at(i);
        if(letter != 'A' && letter != 'C' &&
           letter != 'G' && letter != 'T'  )
        {
            //shc_log_error("comp %u, %u read contain unknown char at %d, %s\n",
            //                curr_comp, i_read, i, read.c_str());
            //exit(0);
            return false;
        }
    }
    if(read.size()<kmer_length+2)
    {
        return false;
    }
    return true;
}

void Sequence_graph_handler::
load_all_paired_read_no_concat(
    std::string read_path_p1, std::string read_path_p2, Kmer_Node_map & kmer_node_map)
{
    Single_read_list & read_list_p1 = coll_read_list.p1_reads;
    Single_read_list & read_list_p2 = coll_read_list.p2_reads;
    std::string temp, read_base_p1, read_base_p2;
    std::ifstream reads_p1_reader(read_path_p1);
    std::ifstream reads_p2_reader(read_path_p2);

    uint64_t i_th_read = 0;
    uint64_t num_valid_pair = 0;
    read_num_t read_index_p1;
    read_num_t read_index_p2;

    Kmer_counter_map temp_kmer_map;
    //Read_subsampler read_subsampler(subsample_factor);

    while(true)
    {
        if((std::getline(reads_p1_reader, temp, '\n')).eof())
            break;
        else
            std::getline(reads_p1_reader, read_base_p1, '\n');

        if((std::getline(reads_p2_reader, temp, '\n')).eof())
            break;
        else
            std::getline(reads_p2_reader, read_base_p2, '\n');

        if(!trim_read(read_base_p1, i_th_read))
            continue;
        if(!trim_read(read_base_p2, i_th_read))
            continue;

        bool is_valid_pair_vd = false;
        vd_t p1_vd_end = traverse_read_is_all_kmer_node_valid(read_base_p1,
                                        kmer_node_map, PAIR_1);
        vd_t p2_vd_front = NULL;
        if(p1_vd_end!=NULL)
        {
            p2_vd_front = traverse_read_is_all_kmer_node_valid(read_base_p2,
                                        kmer_node_map, PAIR_2);
            if(p2_vd_front!=NULL)
            {
                is_valid_pair_vd = true;
            }
        }
        // add reads

        read_list_p1.add_read(read_base_p1, read_index_p1); //read index local to its dict
        read_list_p2.add_read(read_base_p2, read_index_p2); //read index local to its dict

        // assign terminal nodes for each pair of read
        if(is_valid_pair_vd)
        {
            //shc_log_info(shc_logname, "%u\n", i_th_read);
            //shc_log_info(shc_logname, "p1 %s\n", read_base_p1.c_str());
            //shc_log_info(shc_logname, "p2 %s\n", read_base_p2.c_str());

            graph[p1_vd_end].term_nodes_info.emplace_back(PAIR_END_NODE, i_th_read);
            graph[p2_vd_front].term_nodes_info.emplace_back(PAIR_FRONT_NODE, i_th_read);
            num_valid_pair++;
            p1_end.push_back(p1_vd_end);
            p2_front.push_back(p2_vd_front);
        }
        else
        {
            p1_end.push_back(NULL);
            p2_front.push_back(NULL);
        }

        // link two read regardless
        //coll_read_list.link_two_reads(read_index_p1, read_index_p2);
        i_th_read++;
    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "finish load pair read\n");
#endif
    reads_p1_reader.close();
    reads_p1_reader.clear();
    reads_p2_reader.close();
    reads_p2_reader.clear();
}

void Sequence_graph_handler::load_all_paired_read(Kmer_Node_map & kmer_node_map)
{
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "start load pair read\n");
#endif
    Single_read_list & read_list_p1 = coll_read_list.p1_reads;
    Single_read_list & read_list_p2 = coll_read_list.p2_reads;
    std::string temp, read_base_p1, read_base_p2;

    int p1_i = 0;
    int p2_i = 0;
    std::string read_path_p1 = read_path_p1_prefix + "_" + std::to_string(p1_i);
    std::string read_path_p2 = read_path_p2_prefix + "_" + std::to_string(p2_i);

    std::ifstream reads_p1_reader(read_path_p1);
    std::ifstream reads_p2_reader(read_path_p2);

    //std::ofstream p1_d(setting.local_files.output_path + "/read1");
    //std::ofstream p2_d(setting.local_files.output_path + "/read2");

    uint64_t i_th_read = 0;
    uint64_t num_valid_pair = 0;
    read_num_t read_index_p1;
    read_num_t read_index_p2;

    Kmer_counter_map temp_kmer_map;
    //Read_subsampler read_subsampler(subsample_factor);

    while(true)
    {
        //update_edge_count_with_read(read_base_p1, kmer_node_map);
        //shc_log_info(shc_logname, "%u: read1 %s\n", i_th_read, read_base_p1.c_str());
        //shc_log_info(shc_logname, "%u: read2 %s\n", i_th_read, read_base_p2.c_str());
        if((std::getline(reads_p1_reader, temp, '\n')).eof())
        {
            reads_p1_reader.close();
            reads_p1_reader.clear();
            read_path_p1 = read_path_p1_prefix + "_" + std::to_string(++p1_i);
            //std::cout << "try " << read_path_p1 << std::endl;
            if(exist_path(read_path_p1))
            {
                reads_p1_reader.open(read_path_p1);
                std::getline(reads_p1_reader, temp, '\n');
                std::getline(reads_p1_reader, read_base_p1, '\n');
                //std::cout << "read_base_p2: " << read_base_p1 << std::endl;
            }
            else
            {
                //std::cout << "read 1 break" << std::endl;
                break; // case when no file left to process, and last reach eof
            }

        }
        else
            std::getline(reads_p1_reader, read_base_p1, '\n');

        //p1_d << temp << std::endl << read_base_p1 << std::endl;

        if((std::getline(reads_p2_reader, temp, '\n')).eof())
        {
            reads_p2_reader.close();
            reads_p2_reader.clear();
            read_path_p2 = read_path_p2_prefix + "_" + std::to_string(++p2_i);
            //std::cout << "try " << read_path_p2 << std::endl;
            if(exist_path(read_path_p2))
            {
                reads_p2_reader.open(read_path_p2);
                std::getline(reads_p2_reader, temp, '\n');
                std::getline(reads_p2_reader, read_base_p2, '\n');
                //std::cout << "read_base_p2: " << read_base_p2 << std::endl;
            }
            else
            {
                //std::cout << "read 2 break" << std::endl;
                break; // case when no file left to process, and last reach eof
            }
        }
        else
            std::getline(reads_p2_reader, read_base_p2, '\n');

        //p2_d << temp << std::endl << read_base_p2 << std::endl;

        if(!trim_read(read_base_p1, i_th_read))
            continue;
        if(!trim_read(read_base_p2, i_th_read))
            continue;
        //read_base_p1.resize(setting.pair_1_read_length);
        //read_base_p2.resize(setting.pair_2_read_length);

        bool is_valid_pair_vd = false;
        vd_t p1_vd_end = traverse_read_is_all_kmer_node_valid(read_base_p1,
                                        kmer_node_map, PAIR_1);
        vd_t p2_vd_front = NULL;
        if(p1_vd_end!=NULL)
        {
            p2_vd_front = traverse_read_is_all_kmer_node_valid(read_base_p2,
                                        kmer_node_map, PAIR_2);
            if(p2_vd_front!=NULL)
            {
                is_valid_pair_vd = true;
            }
        }
        // add reads

        //if(read_subsampler.simple_decider())
        //{
            read_list_p1.add_read(read_base_p1, read_index_p1); //read index local to its dict
            read_list_p2.add_read(read_base_p2, read_index_p2); //read index local to its dict

            // assign terminal nodes for each pair of read
            if(is_valid_pair_vd)
            {
                //shc_log_info(shc_logname, "%u\n", i_th_read);
                //shc_log_info(shc_logname, "p1 %s\n", read_base_p1.c_str());
                //shc_log_info(shc_logname, "p2 %s\n", read_base_p2.c_str());

                graph[p1_vd_end].term_nodes_info.emplace_back(PAIR_END_NODE, i_th_read);
                graph[p2_vd_front].term_nodes_info.emplace_back(PAIR_FRONT_NODE, i_th_read);
                num_valid_pair++;
                p1_end.push_back(p1_vd_end);
                p2_front.push_back(p2_vd_front);
            }
            else
            {
                p1_end.push_back(NULL);
                p2_front.push_back(NULL);
            }
        //}

        // link two read regardless
        //coll_read_list.link_two_reads(read_index_p1, read_index_p2);
        i_th_read++;
    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "finish load pair read\n");
#endif
    reads_p1_reader.close();
    reads_p1_reader.clear();
    reads_p2_reader.close();
    reads_p2_reader.clear();
}


vd_t Sequence_graph_handler::
traverse_read_is_all_kmer_node_valid(std::string & base,
                           Kmer_Node_map & kmer_node_map, uint8_t read_type)
{
    Kmer_Node_map_iterator it_end = kmer_node_map.end();
    std::string begin_kmer(base.substr(0, kmer_length));
    std::string last_kmer(base.substr(base.size()-kmer_length, kmer_length));

    Kmer_Node_map_iterator it_begin = kmer_node_map.find(begin_kmer);
    Kmer_Node_map_iterator it_curr = it_begin;
    Kmer_Node_map_iterator it_last = kmer_node_map.find(last_kmer);
    if(it_curr==it_end)
    {
        //shc_log_info(shc_logname, "no beginning kmer %s\n", begin_kmer.c_str());
        return NULL;
    }
    if( it_last==it_end)
    {
        //shc_log_info(shc_logname, "no last kmer %s\n", last_kmer.c_str());
        return NULL;
    }
    for(int i=1; i<base.size()-kmer_length; i++)
    {
        std::string kmer(base.substr(i, kmer_length));
        it_curr = kmer_node_map.find(kmer);
        if(it_curr == it_end)
        {
            //shc_log_info(shc_logname, "cannot find %s\n", kmer.c_str());
            return NULL;
        }
    }
    if(read_type==PAIR_1)
        return it_last->second;
    if(read_type==PAIR_2)
        return it_begin->second;
}

void Sequence_graph_handler::
load_all_single_read(std::string& read_path, Kmer_Node_map & kmer_node_map)
{
    Single_read_list & read_list = coll_read_list.s_reads;
    std::string kmer_base, temp, read_base;
    std::ifstream read_file_reader(read_path);
    kmer_base.resize(kmer_length);

    int i = 0;
    int j = 0;
    read_num_t read_index;

    Kmer_counter_map temp_kmer_map;
    //Read_subsampler read_subsampler(subsample_factor);

    while(  std::getline(read_file_reader, temp) &&
            std::getline(read_file_reader, read_base))
    {
        //std::cout << "temp " << temp << std::endl;
        //std::cout << "read_base " << read_base << std::endl;
        //shc_log_info(shc_logname, "single reads %s\n", read_base.c_str());
        //update_edge_count_with_read(read_base, kmer_node_map);
        if(!trim_read(read_base, i))
            continue;
        //if(read_subsampler.simple_decider())
        //{
            read_list.add_read(read_base, read_index);
            j++;
        //}
        i++;
    }
    read_file_reader.close();
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "load %d single reads\n", j);
#endif
}

void Sequence_graph_handler::
use_pair_reads_to_build_edge(std::string& read_path_1, std::string& read_path_2,
     Kmer_Node_map& kmer_node_map, read_num_t *read_id_ptr)
{
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start using pair reads to build graph\n");
#endif
    Single_read_list & read_list_p1 = coll_read_list.p1_reads;
    Single_read_list & read_list_p2 = coll_read_list.p2_reads;
    std::string temp, read_base_p1, read_base_p2;
    std::ifstream reads_p1_reader(read_path_1);
    std::ifstream reads_p2_reader(read_path_2);

    while(  std::getline(reads_p1_reader, temp) &&
            std::getline(reads_p1_reader, read_base_p1) &&
            std::getline(reads_p2_reader, temp) &&
            std::getline(reads_p2_reader, read_base_p2))
    {
        vd_t vd_end = process_one_read_to_build_graph(
                kmer_node_map, read_base_p1, PAIR_END_NODE);
        vd_t vd_front = process_one_read_to_build_graph(
                kmer_node_map, read_base_p2, PAIR_FRONT_NODE);
        read_num_t read_index1, read_index2;
        if(read_list_p1.add_read(read_base_p1, read_index1))
            (*read_id_ptr)++;
        if(read_list_p2.add_read(read_base_p2, read_index2))
            (*read_id_ptr)++;
        //coll_read_list.link_two_reads(read_index1, read_index2);

        if(vd_end!=NULL && vd_front!=NULL)
        {
            //paired_node_set.insert(vd_end);
            //paired_node_set.insert(vd_front);
            p1_end.push_back(vd_end);
            p2_front.push_back(vd_front);
            graph[vd_end].term_nodes_info.emplace_back(PAIR_END_NODE, read_index1);
            graph[vd_front].term_nodes_info.emplace_back(PAIR_FRONT_NODE, read_index2);
        }
        else
        {
            p1_end.push_back(NULL);
            p2_front.push_back(NULL);
        }
    }
    reads_p1_reader.close();
    reads_p2_reader.close();
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start using pair reads to build graph\n");
#endif
}


void Sequence_graph_handler::
use_reads_to_build_edge(std::string& read_path,
                        Kmer_Node_map& kmer_node_map, read_num_t *read_id_ptr)
{
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start using single read to build edge\n");
#endif
    Single_read_list & read_list = coll_read_list.s_reads;
    std::string kmer_base, temp, read_base;
    std::ifstream read_file_reader(read_path);
    kmer_base.resize(kmer_length);

    while(  std::getline(read_file_reader, temp) &&
            std::getline(read_file_reader, read_base))
    {
        process_one_read_to_build_graph( kmer_node_map, read_base, NO_PAIR);
        read_num_t read_index;
        if(read_list.add_read(read_base, read_index))
            (*read_id_ptr)++;
    }
    read_file_reader.close();
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Finish using single read to build edge\n");
#endif
}

void Sequence_graph_handler::assign_read_to_xnode()
{
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start assign reads to x node\n");
#endif
    start_timer(&timer);
    typedef tsl::hopscotch_map<uint64_t, std::vector<vd_t>,
                   hash_u64, equ64> Kmer_vds_map;
    typedef tsl::hopscotch_map<uint64_t, std::vector<vd_t>,
                    hash_u64, equ64>::iterator Kmer_vds_map_iterator;
    Kmer_vds_map kmer_vds_map;
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "constructing x-node k start seq map\n");
#endif
    //for(std::set<vd_t>::iterator it=xnode_set.begin(); it!=xnode_set.end(); it++)
    //{
    //    vd_t vd = *it;
    //    std::string k_start(graph[vd].seq.substr(0,kmer_length));
    //    start_seq_node_mmap.insert(std::make_pair(k_start, vd));
    //}

    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    uint64_t byte;
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        encode_kmer(&(graph[vd].seq.at(0)), &byte, kmer_length);
        Kmer_vds_map_iterator kv_it = kmer_vds_map.find(byte);
        if(kv_it != kmer_vds_map.end())
        {
            (kv_it.value()).push_back(vd);
        }
        else
        {
            std::vector<vd_t> vec_seq(1, vd);
            //vec_seq.reserve(KMER_VEC_INIT_SIZE);
            kmer_vds_map.insert(std::make_pair(byte, vec_seq));
        }
    }

    //examine each subseq
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Reading lines\n");
#endif
    num_assigned_read = 0;
    //std::ofstream debug_read_file(setting.local_files.output_path+"/all_read");
    for(size_t j = 0; j < coll_read_list.get_num_reads(); j++)
    {
        Read_acc acc = coll_read_list.get_read(j);
        //shc_log_info(shc_logname, "read len %u\n", acc.len);
        //std::string a_read(acc.read_ptr, acc.len);
        //shc_log_info(shc_logname, "%s has count %u\n", a_read.c_str(), acc.read_count);
        if(acc.len > kmer_length)
        {
            for(int i=1; i<acc.len - kmer_length; i++)
            {
                encode_kmer(acc.read_ptr + i, &byte, kmer_length);
                Kmer_vds_map_iterator it = kmer_vds_map.find(byte);
                if(it != kmer_vds_map.end())
                {
                    std::vector<vd_t> & vd_vec = it.value();
                    for(std::vector<vd_t>::iterator vd_it =vd_vec.begin();
                                                    vd_it!=vd_vec.end(); ++vd_it)
                    {
                        if(is_read_bridge_node(acc.read_ptr, acc.len, i, *vd_it))
                        {
                            num_assigned_read ++;
                            node_link_read(*vd_it, j, i); //j is virtual read index
                        }
                    }
                }
            }
        }
    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Finish assign reads to x node\n");
#endif
#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "load and match reads finish, ";
    stop_timer(&timer);
#endif
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
    return memcmp(&(graph[vd].seq.at(0)),
            read_ptr+info_start, graph[vd].seq_len())==0;
}


uint64_t Sequence_graph_handler::
resolve_all_pair_reads(int search_hop, std::set<vd_t> & check_vd)
{
    assert(p1_end.size()==p2_front.size());
    //std::cout << "Start resolve all pairs\n" << std::endl;
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start resolve all pairs\n");
#endif
    uint64_t num_pair_read = p1_end.size();

    Node_pair_path_map node_pair_path_map;
    std::set<Node_pair> checked_read_pair;

    uint64_t num_useful_sr = 0;
    uint64_t num_sr = 0;
    for(uint64_t i=0; i<num_pair_read; i++)
    {
        //shc_log_info(shc_logname, "\n");

#ifdef LOG_PAIR_SEARCH
        shc_log_info(shc_logname, "\t\t comp %u align read index %u\n", curr_comp, i);
#endif
        vd_t vd1 = p1_end[i];
        vd_t vd2 = p2_front[i];


        if(vd1!=NULL && vd2!=NULL)
        {

            bool is_sr = false;
            if(check_and_add_simple_read_pair(i, is_sr))
            {
                num_useful_sr++;
                check_vd.insert(vd1); // vd1 == vd2
            }

            if(is_sr)
                num_sr++;


            if(search_hop > 0 && !is_sr) // reads are not simply connected
            {
                // pair has not been checked yet
                if(checked_read_pair.find(Node_pair(vd1, vd2))==checked_read_pair.end())
                {
                    std::vector<vd_t> path;

                    read_num_t l_start_read_id =
                            coll_read_list.convert_align_index_to_local_index(i, true);

                    if(search_path_in_pair_nodes_with_bfs(
                            vd1, vd2, search_hop, path, node_pair_path_map, l_start_read_id))
                    {
#ifdef LOG_PAIR_SEARCH
                        info_log_info(shc_logname, "comp %d path ", curr_comp);
                        for(std::vector<vd_t>::iterator path_it= path.begin();
                                               path_it!=path.end(); path_it++)
                        {
                            info_log_info(shc_logname, "%u ", graph[*path_it].node_id);
                            check_vd.insert(*path_it);
                        }
                        info_log_info(shc_logname, "\n");
#endif
                        std::string middle;

                        read_num_t l_stop_read_id =
                                coll_read_list.convert_align_index_to_local_index(i, false);
                        read_num_t v_start_read_id =
                                coll_read_list.convert_local_index_to_virtual_index(l_start_read_id, true);
                        read_num_t v_stop_read_id =
                                coll_read_list.convert_local_index_to_virtual_index(l_stop_read_id, false);

                        convert_vd_path_to_string(path, middle, l_start_read_id, l_stop_read_id);
                        clear_term_info(i);
                        if(middle.size() > setting.seq_graph_setup.mate_pair_len)
                        {
#ifdef LOG_PAIR_SEARCH
                            shc_log_info(shc_logname,  "middle too long %d\n", middle.size());
#endif
                            continue;
                        }
                        else
                        {
                            shc_log_info(shc_logname,  "good non-simple super read %d\n", middle.size());
                        }
#ifdef LOG_PAIR_SEARCH
                        shc_log_info(shc_logname, "comp %d middle is %s\n", curr_comp, middle.c_str());
#endif
                        coll_read_list.add_super_read(v_start_read_id, v_stop_read_id,
                                          l_start_read_id, l_stop_read_id, middle);

                        num_sr++;
                        num_useful_sr++;
#ifdef LOG_PAIR_SEARCH
                        Super_read sr(l_start_read_id, l_stop_read_id, middle);
                        Read_acc acc = coll_read_list.get_super_read(sr);
                        std::string sr_str(acc.read_ptr, acc.len);
                        //shc_log_info(shc_logname, "SR:\n");
                        //log_read_seq(acc);

                        //read and terminal node
                        bundled_node_p & p1_node = graph[path.front()];
                        bundled_node_p & p2_node = graph[path.back()];
                        Read_acc p1_acc = coll_read_list.get_read(v_start_read_id);
                        Read_acc p2_acc = coll_read_list.get_read(v_stop_read_id);
                        std::string read1(p1_acc.read_ptr, p1_acc.len);
                        std::string read2(p2_acc.read_ptr, p2_acc.len);

                        std::size_t p1_found = read1.find(p1_node.seq);
                        if (p1_found != std::string::npos)
                        {
                            node_link_read(path.front(), v_start_read_id, p1_found); //j is virtual read index
                            //shc_log_info(shc_logname, "left term node included by the left read\n");
                        }

                        std::size_t p2_found = read2.find(p2_node.seq);
                        if (p2_found != std::string::npos)
                        {
                            node_link_read(path.back(), v_start_read_id, p2_found); //j is virtual read index
                            //shc_log_info(shc_logname, "right term node included by the right read\n");
                        }

                        for(std::vector<vd_t>::iterator path_it= path.begin()+1;
                                               path_it!=path.end()-1; path_it++)
                        {
                            bundled_node_p path_node = graph[*path_it];
                            //shc_log_info(shc_logname, "Node:\n");
                            //shc_log_info(shc_logname, "%s\n", path_node.seq.c_str());

                            std::size_t found = sr_str.find(path_node.seq);
                            if (found!=std::string::npos)
                            {
                                //shc_log_info(shc_logname, "at index %u\n", found);
                                node_link_read(*path_it, v_start_read_id, found); //j is virtual read index
                            }
                            else
                            {
                                shc_log_error("path node not covered by SR\n");
                                exit(1);
                            }
                        }

                        shc_log_info(shc_logname, "%u, find non-simple sr %s\n", curr_comp, sr_str.c_str());
#endif
                    }
#ifdef LOG_PAIR_SEARCH
                    else
                        shc_log_info(shc_logname, "%u, does not find a path\n", curr_comp);
#endif
                    checked_read_pair.insert(Node_pair(vd1, vd2));
                }
#ifdef LOG_PAIR_SEARCH
                else
                {
                    shc_log_info(shc_logname, "already processed non-simple pair\n");
                }
#endif
            }

        }
        //shc_log_info(shc_logname, "\n");
    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Finish resolve %u pairs, %u is useful\n",
                                        num_sr, num_useful_sr);
#endif
    return num_useful_sr;
}

/**
 * The function paints the gray node to black when vd check is ancestrial of
 * some node in the queue. Black nodes in the queue will not be processed
 * by BFS
 *
 * Function deals with the case when there are not unique path, but the second
 * path is found very late, and there are already many nodes in the queue
 *  x -> x -> x -> x
 *  |              |
 *  v              v
 *  x <- x <- x <- x
 *  |
 *  X
 * @param vd_root is the vd starting point
 * @param vd_check       vd being compared
 * @param vd_queue
 */
void Sequence_graph_handler::
paint_queue_node(vd_t vd_root, vd_t vd_check, std::deque<vd_t> & vd_queue)
{
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Strat paint nodes in queue\n");
#endif
    for(std::deque<vd_t>::iterator it=vd_queue.begin(); it!=vd_queue.end(); ++it)
    {
        if(is_ancestrial(vd_root, vd_check, *it))
        {
#ifdef LOG_SEQ_GRAPH
            shc_log_info(shc_logname, "%u is painted\n", graph[*it].node_id);
#endif
            graph[*it].s_info.color = BLACK;
            graph[*it].s_info.parent = NULL;
            graph[*it].s_info.d = 0;
            graph[*it].s_info.mid_seq_len = -1;
        }
#ifdef LOG_SEQ_GRAPH
        else
        {
            shc_log_info(shc_logname, "%u is not painted\n", graph[*it].node_id);
        }
#endif
    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Finish paint nodes in queue\n");
#endif
}

// see if vd_check is ancestrial of vd
bool Sequence_graph_handler::
is_ancestrial(vd_t vd_root, vd_t vd_check, vd_t vd)
{
    vd_t parent = static_cast<vd_t>(graph[vd].s_info.parent);
    while(parent != vd_root)
    {
        if(parent==NULL)
            return false;

        if(parent == vd_check)
            return true;
        else
            parent = static_cast<vd_t>(graph[parent].s_info.parent);
    }
    return false;
}

void Sequence_graph_handler::
create_path_vd_set(vd_t vd_end, vd_t vd_start, std::set<vd_t> & vd_set)
{
    vd_set.insert(vd_end);
    vd_t prev_vd = vd_end;
    vd_t parent = static_cast<vd_t>(graph[vd_end].s_info.parent);
    if(parent==NULL)
    {
        shc_log_error("find a NULL parent\n");
        shc_log_info(shc_logname, "%u, %u has parent NULL\n", curr_comp, graph[prev_vd].node_id);
        exit(1);
    }
#ifdef LOG_SEQ_GRAPH
    else
    {
        shc_log_info(shc_logname, "%u, %u has parent %u\n", curr_comp,
                                graph[prev_vd].node_id, graph[parent].node_id);
    }
#endif

    while(parent != vd_start)
    {
        vd_set.insert(parent);
        prev_vd = parent;
        parent = static_cast<vd_t>(graph[parent].s_info.parent);
        if(parent==NULL)
        {
            shc_log_error("find a NULL parent\n");
            shc_log_info(shc_logname, "%u, %u has parent NULL\n", curr_comp, graph[prev_vd].node_id);
            exit(1);
        }
#ifdef LOG_SEQ_GRAPH
        else
        {
            shc_log_info(shc_logname, "%u, %u has parent %u\n", curr_comp,
                                    graph[prev_vd].node_id, graph[parent].node_id);
        }
#endif
    }
}

bool Sequence_graph_handler::
check_and_create_path_vector(vd_t vd_curr, vd_t vd_start, std::vector<vd_t> & vd_vec)
{
    vd_t repeat_vd;
    if(!is_unique_path_to_u(vd_curr, repeat_vd))
    {
        //shc_log_info(shc_logname, "is not unique path to %u\n", graph[vd_curr].node_id);
        return false;
    }

    vd_t parent = static_cast<vd_t>(graph[vd_curr].s_info.parent);
    if (vd_curr == vd_start)
    {
        vd_vec.push_back(vd_curr);
        //shc_log_info(shc_logname, "find a path\n");
        return true;
    }
    else if (parent == NULL)
    {
        //shc_log_info(shc_logname, "No path between two nodes, because a node in the middle is later found to have two path\n");
        return false;
    }
    else
    {
        if(check_and_create_path_vector(parent, vd_start, vd_vec))
        {
            vd_vec.push_back(vd_curr);
            return true;
        }
        else
        {
            return false;
        }
    }
}

void Sequence_graph_handler::
convert_vd_path_to_string(std::vector<vd_t> & path, std::string & middle,
                            read_num_t local_read_1, read_num_t local_read_2)
{
    if(path.size()==0 || path.size()==1)
    {
        shc_log_error("comp %d, vd path contains 0 or 1 node\n", curr_comp);
        exit(1);
    }
    middle.clear();
    vd_t vd_1 = path.at(0);
    vd_t vd_2 = path.back();
    std::string & vd_1_seq = graph[vd_1].seq;
    std::string & vd_2_seq = graph[vd_2].seq;

    Read_acc p1_acc = coll_read_list.p1_reads.get_read(local_read_1);
    std::string read_p1_last(p1_acc.read_ptr+p1_acc.len-kmer_length, kmer_length);
    Read_acc p2_acc = coll_read_list.p2_reads.get_read(local_read_2);
    std::string read_p2_first(p2_acc.read_ptr, kmer_length);

    //shc_log_info(shc_logname, "node %u  seq %s\n", curr_node.node_id, seq.c_str());
    //shc_log_info(shc_logname, "p1 last  seq %s\n", read_p1_last.c_str());
    //shc_log_info(shc_logname, "p2 first seq %s\n", read_p2_first.c_str());
    //shc_log_info(shc_logname, "node 1 seq %s\n", vd_1_seq.c_str());
    //shc_log_info(shc_logname, "node 2 seq %s\n", vd_2_seq.c_str());

    int64_t p1_end_index  = find_terminal_index( vd_1,
            vd_1_seq, p1_acc.read_ptr+p1_acc.len-kmer_length, kmer_length, p1_acc.len, true);

    int64_t p2_first_index =  find_terminal_index(vd_2,
                     vd_2_seq, p2_acc.read_ptr, kmer_length, p2_acc.len, false);

    std::string first_vd_last = vd_1_seq.substr(p1_end_index);
    std::string last_vd_begin = vd_2_seq.substr(0, p2_first_index);
    middle += first_vd_last;
    //shc_log_info(shc_logname, "middle %s\n", middle.c_str());
    for(std::vector<vd_t>::iterator it=path.begin()+1; it!=path.end()-1; ++it)
    {
        bundled_node_p & node = graph[*it];
        aer_t aer = boost::edge(*(it-1),*(it), graph);
        if(!aer.second)
        {
            shc_log_error("edge in sr path does not exist\n");
            exit(0);
        }
        bundled_edge_p & edge_p = graph[aer.first];
        //shc_log_info(shc_logname, "edge weight %d\n", edge_p.weight);
        //shc_log_info(shc_logname, "node seq %s\n", node.seq.c_str());

        middle += node.seq.substr(edge_p.weight);
        //shc_log_info(shc_logname, "middle %s\n", middle.c_str());
    }

    vd_t vd_second_last = path.at(path.size()-2);
    aer_t aer = boost::edge(vd_second_last, vd_2, graph);
    if(!aer.second)
    {
        shc_log_error("edge in sr path does not exist\n");
        exit(0);
    }

    //shc_log_info(shc_logname, "last_vd_begin %s\n", last_vd_begin.c_str());
    bundled_edge_p & edge_p = graph[aer.first];
    //shc_log_info(shc_logname, "last edge weight %d\n", edge_p.weight);
    if(middle.size() > edge_p.weight)
    {
        middle.resize(middle.size() - edge_p.weight);
        middle += last_vd_begin;
    }
    else
    {
        int64_t index = std::min(edge_p.weight, last_vd_begin.size());
        middle += last_vd_begin.substr(0, index);
        //shc_log_error("check for strange case\n");
        //exit(1);
    }
}

bool Sequence_graph_handler::is_unique_path_to_u(vd_t u, vd_t & repeat_vd)
{
    in_eip_t in_eip = boost::in_edges(u, graph);
    int num_black_num = 0;
    for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
    {
        vd_t vd_source = boost::source(*in_ei, graph);
        if(graph[vd_source].s_info.color == GRAY)
        {
            repeat_vd = u;
            //shc_log_info(shc_logname, "a gray node pointed by another gray node\n");
            return false;
        }

        // pointed by more than one black nodes
        if(graph[vd_source].s_info.color == BLACK)
        {
            if(num_black_num == 0)
            {
                num_black_num++;
            }
            else
            {
                //shc_log_info(shc_logname, "pointed by %d black nodes\n", num_black_num);
                repeat_vd = u;
                //shc_log_info(shc_logname, "repeat node %u\n", graph[repeat_vd].node_id);
                return false;
            }
            // or the there is a loop between two nodes
            aer_t aer = boost::edge(u, vd_source, graph);
            if(aer.second) // contain such edge
            {
                //shc_log_info(shc_logname, "contain two node loop, %u %u\n",
                //               graph[u].node_id, graph[vd_source].node_id);
                repeat_vd = vd_source;
                //shc_log_info(shc_logname, "repeat node %u\n", graph[repeat_vd].node_id);
                return false;
            }
        }
    }
    return true;
}

void Sequence_graph_handler::reset_node_color(std::set<vd_t> & colored_set)
{
    for(std::set<vd_t>::iterator it = colored_set.begin();
                                 it!= colored_set.end(); it++)
    {
        graph[*it].s_info.color = WHITE;
        graph[*it].s_info.parent = NULL;
        graph[*it].s_info.mid_seq_len = 0;
        graph[*it].s_info.d = 0;
    }
    colored_set.clear();
}

void Sequence_graph_handler::
get_all_descendant(std::set<vd_t> & repeat_set, vd_t root_vd)
{
    //shc_log_info(shc_logname, "start get all descnednat, size %d\n", repeat_set.size());
    std::map<vd_t, uint8_t> color_map;
    std::deque<vd_t> gray_queue;
    gray_queue.push_back(root_vd);
    color_map[root_vd] = GRAY;
    while (gray_queue.size()!= 0)
    {
        vd_t u = gray_queue.front();
        //shc_log_info(shc_logname, "gray queue pushed %d\n", graph[u].node_id);
        repeat_set.insert(u);
        gray_queue.pop_front();
        out_eip_t out_eip = boost::out_edges(u, graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t v = boost::target(*out_ei, graph);
            if(color_map.find(v)==color_map.end())
            {
                color_map[v] = GRAY;
                gray_queue.push_back(v);
            }
        }
    }
    color_map.clear();
    //shc_log_info(shc_logname, "finish get all descnednatsize %d\n", repeat_set.size());
}




/**
 * Function uses BFS from vd1 to vd2
 * @param vd1
 * @param vd2
 * @param curr_hop
 * @return True if there is a unique path after exhausting search, path store
 *              inside argument path, (inclusive on both ends)
 *         False 1. if no path found within max hop distance or
 *               2. if the path is no unique
 */
bool Sequence_graph_handler::
search_path_in_pair_nodes_with_bfs(vd_t vd1, vd_t vd2, int max_hop,
                   std::vector<vd_t> & path,
                   Node_pair_path_map & node_pair_path_map, read_num_t local_read_1)
{
#ifdef LOG_PAIR_SEARCH
    shc_log_info(shc_logname, "%u, Start read pair search for %u %u\n",
                                curr_comp, graph[vd1].node_id, graph[vd2].node_id);
#endif

    std::set<vd_t> colored_set;
    std::set<vd_t> repeat_vds;


    std::deque<vd_t> vd_queue;
    std::set<vd_t> st_path_vd_set; // source target path vd set
    bool is_find_first = false;
    bool is_find_second = false;
    graph[vd1].s_info.color = GRAY;

    // find mid seq len
    std::string & vd_1_seq = graph[vd1].seq;
    Read_acc p1_acc = coll_read_list.p1_reads.get_read(local_read_1);
    std::string read_p1_last(p1_acc.read_ptr+p1_acc.len-kmer_length, kmer_length);
    int64_t p1_end_index  = find_terminal_index( vd1,
            vd_1_seq, p1_acc.read_ptr+p1_acc.len-kmer_length, kmer_length, p1_acc.len, true);
    read_length_t mid_seq_len = graph[vd1].seq_len() - p1_end_index;
    shc_log_info(shc_logname, "init mid_seq_len %d\n", mid_seq_len);
    graph[vd1].s_info.mid_seq_len = mid_seq_len;


    vd_queue.push_back(vd1);

    while (vd_queue.size()>0)
    {
        vd_t u = vd_queue.front();
        colored_set.insert(u);

        vd_queue.pop_front();
        S_info & u_bfs_info = graph[u].s_info;

        if(u_bfs_info.mid_seq_len > setting.seq_graph_setup.mate_pair_len)
        {
            shc_log_info(shc_logname, "too long %d\n", u_bfs_info.mid_seq_len);

            continue;
        }

#ifdef LOG_PAIR_SEARCH
        shc_log_info(shc_logname, "%u, curr element is %d\n", curr_comp, graph[u].node_id);
        shc_log_info(shc_logname, "%u, after pop queue ", curr_comp);
        for(std::deque<vd_t>::iterator queue_it= vd_queue.begin();
                                       queue_it!=vd_queue.end(); queue_it++)
        {
            info_log_info(shc_logname, "%u ", graph[*queue_it].node_id);
        }
        info_log_info(shc_logname, "\n");
#endif
        //Fail condition 1
        if(u_bfs_info.d > max_hop && !is_find_first)
        {
#ifdef LOG_PAIR_SEARCH
            shc_log_info(shc_logname,
              "%u, No resolve,searched greater than max hop, but do not find a path\n", curr_comp);
#endif
            reset_node_color(colored_set);
            return false;
        }
        if (u_bfs_info.color == WHITE)
        {
            shc_log_error("comp %d, impossible white nodes", curr_comp);
            exit(1);
        }
        else if(u_bfs_info.color == BLACK)
        {
            // when some nodes are painted
#ifdef LOG_PAIR_SEARCH
            shc_log_info(shc_logname, "%u, a black node in the queue, skip\n", curr_comp);
#endif
            continue;
        }
        // when this gray nodes is pointed by other gray nodes, or by two black nodes
        vd_t vd_repeat;
        bool is_unique_path = is_unique_path_to_u(u, vd_repeat);
        if(!is_unique_path)
        {
            u_bfs_info.color = BLACK;
            repeat_vds.insert(vd_repeat);
            get_all_descendant(repeat_vds, vd_repeat);
            if(repeat_vds.find(vd2) != repeat_vds.end())
            {
#ifdef LOG_PAIR_SEARCH
                shc_log_info(shc_logname, "%u, target inside repeat node\n", curr_comp);
#endif
                return false;
            }
#ifdef LOG_PAIR_SEARCH
            shc_log_info(shc_logname, "%u, find non-unique path, node %u is colored to be black\n",
                               curr_comp, graph[u].node_id);
#endif
            paint_queue_node(vd1, vd_repeat, vd_queue);
            if(is_find_first)
            {
                if(is_ancestrial(vd1, vd_repeat, vd2))
                {
#ifdef LOG_PAIR_SEARCH
                    shc_log_info(shc_logname, "target node visit a node which has cycle\n");
#endif
                    return false;
                }
            }

        }
        else
        {
#ifdef LOG_PAIR_SEARCH
            shc_log_info(shc_logname, "%u, find a unique path to %u\n", curr_comp, graph[u].node_id);
#endif
            // paint gray to out edges nodes from u
            out_eip_t out_eip = boost::out_edges(u, graph);
            for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
            {
                vd_t v = boost::target(*out_ei, graph);
#ifdef LOG_PAIR_SEARCH
                shc_log_info(shc_logname, "%u, explore v %u\n", curr_comp, graph[v].node_id);
#endif
                if(st_path_vd_set.find(v) != st_path_vd_set.end())
                {
#ifdef LOG_PAIR_SEARCH
                    shc_log_info(shc_logname, "%u, not unqiue in the middle, false\n", curr_comp);
#endif
                    reset_node_color(colored_set);
                    return false;
                }

                if(repeat_vds.find(v)!=repeat_vds.end())
                {
#ifdef LOG_PAIR_SEARCH
                    shc_log_info(shc_logname, "%u, on gray node, see a repeat node\n", curr_comp);
#endif
                    continue;
                }



                S_info & v_bfs_info = graph[v].s_info;
                if(v == vd2) // find a path to target
                {
                    if(is_find_first)
                    {
                        is_find_second = true;
#ifdef LOG_PAIR_SEARCH
                        shc_log_info(shc_logname, "not unqiue at the end, false\n");
#endif
                        reset_node_color(colored_set);
                        return false;   // verified path not unique
                    }
                    else
                    {
#ifdef LOG_PAIR_SEARCH
                        shc_log_info(shc_logname, "comp %d, first time find target\n", curr_comp);
                        shc_log_info(shc_logname, "%u has parent %u\n",
                                graph[v].node_id, graph[u].node_id);
#endif
                        v_bfs_info.parent = (void*)u;
                        v_bfs_info.color = WHITE;
                        v_bfs_info.d++;
                        create_path_vd_set(vd2, vd1, st_path_vd_set);
                        is_find_first = true;
                        continue;
                    }
                }

                if (v_bfs_info.color == WHITE) //if is white, regular
                {
#ifdef LOG_PAIR_SEARCH
                    shc_log_info(shc_logname, "gray node %u points to white node %u\n",
                                graph[u].node_id, graph[v].node_id);
                    shc_log_info(shc_logname, "%u has parent %u\n",
                                graph[v].node_id, graph[u].node_id);
#endif
                    v_bfs_info.color = GRAY;
                    v_bfs_info.parent = (void*)u;
                    v_bfs_info.d++;
                    v_bfs_info.mid_seq_len += graph[v].seq_len()-graph[*out_ei].weight + graph[u].s_info.mid_seq_len;
#ifdef LOG_PAIR_SEARCH
                    shc_log_info(shc_logname, "mid seq len %d\n", v_bfs_info.mid_seq_len);
#endif
                    vd_queue.push_back(v);
                }
                else if(v_bfs_info.color == GRAY)
                {
#ifdef LOG_PAIR_SEARCH
                    shc_log_info(shc_logname, "gray node %u points to other gray node %u\n",
                                     graph[u].node_id, graph[v].node_id);
#endif
                    v_bfs_info.color = BLACK;
                    break;
                }
                else if(v_bfs_info.color == BLACK)
                {
                    // a pretty bad case
#ifdef LOG_PAIR_SEARCH
                    shc_log_info(shc_logname, "gray node %u points to black node %u\n",
                                       graph[u].node_id, graph[v].node_id);
                    shc_log_info(shc_logname, "%u has parent NULL\n",
                                graph[v].node_id);
#endif
                    repeat_vds.insert(v);
                    get_all_descendant(repeat_vds, v);
                    if(repeat_vds.find(vd2) != repeat_vds.end())
                    {
#ifdef LOG_PAIR_SEARCH
                        shc_log_info(shc_logname, "%u, target inside repeat node\n", curr_comp);
#endif
                        return false;
                    }
                    if(is_ancestrial(vd1, v, u))
                    {
#ifdef LOG_PAIR_SEARCH
                        shc_log_info(shc_logname, "%u is ancestrial of %u\n",
                                       graph[v].node_id, graph[u].node_id);
#endif
                        paint_queue_node(vd1, v, vd_queue);
                        v_bfs_info.parent = (void*)(NULL);
                        break;
                    }

                    paint_queue_node(vd1, v, vd_queue);
                    v_bfs_info.parent = (void*)(NULL);


                }
                else
                {
                    shc_log_error("comp %d, unknown color\n", curr_comp);
                    exit(1);
                }
            }
            u_bfs_info.color = BLACK;
#ifdef LOG_PAIR_SEARCH
            if(exist_vd.find(u)==exist_vd.end())
                shc_log_info(shc_logname, "vertex u do not exist in exist_vd\n");

            shc_log_info(shc_logname, "%u, gray node is colored to be black %u\n",
                                        curr_comp, graph[u].node_id);
#endif
        }
    }

    if(!is_find_second && is_find_first)
    {
#ifdef LOG_PAIR_SEARCH
        shc_log_info(shc_logname, "pas the first test, wait to back-check\n");
#endif
        if(check_and_create_path_vector(vd2, vd1, path))
        {
            node_pair_path_map.insert(std::make_pair(Node_pair(vd1,vd2), path));
            reset_node_color(colored_set);
#ifdef LOG_PAIR_SEARCH
            shc_log_info(shc_logname, "useful non-simple sr\n");
#endif
            return true;
        }
        reset_node_color(colored_set);
        return false;
    }
    else
    {
#ifdef LOG_PAIR_SEARCH
        shc_log_info(shc_logname, "No path from vd1 to vd2 due to no element in"
                "queu, and search does not reach to specified %d depth\n", max_hop);
#endif
        reset_node_color(colored_set);
        return false;
    }
}

/**
 * This function creates a new node which capture vd1->vd2
 * @param vd1
 * @param vd2
 */
/*
vd_t Sequence_graph_handler::
merge_two_vertices(vd_t vd1, vd_t vd2)
{
    shc_log_info(shc_logname, "merging two vertices\n");
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

    double merged_count = (node1.count + node2.count)/2;
    vd_t vd = boost::add_vertex(graph);
    graph[vd] = bundled_node_p(merged_seq, merged_count, glob_node_id++);

    //merge reads, if the reads are not empty
    if(!node1.reads_info.empty())
        graph[vd].reads_info = node1.reads_info;

    shc_log_info(shc_logname, "before link nodes\n");
    log_node_info(vd, true , false, false);
    log_node_info(vd1, true , false, false);
    //shc_log_info(shc_logname, "update links\n");
    link_in_edges(vd1, vd);
    link_out_edges(vd2, vd);
    shc_log_info(shc_logname, "after link nodes\n");
    log_node_info(vd, true , false, false);
    log_node_info(vd1, true , false, false);


    //update_terminal_node(vd1, vd);
    //update_terminal_node(vd2, vd);


    shc_log_info(shc_logname, "Finish merge vertices\n");
    return vd;
}
*/


/**
*  This function update read terminal vectors. Since the old terminal node (vd1)
*  is replaced with the new merged node (vd). In addition, it copies the
*  terminal node information from old node to new node.
*/
void Sequence_graph_handler::
transfer_read_align_index(vd_t vd_old, vd_t vd_new)
{

    std::vector<Terminal_node_info> & v1_term_nodes_info =
                                            graph[vd_old].term_nodes_info;
    //assert(!v1_term_nodes_info.empty());
    for(std::vector<Terminal_node_info>::iterator
        it =  v1_term_nodes_info.begin();
        it != v1_term_nodes_info.end(); it++)
    {
        assert(PAIR_END_NODE == it->node_type || it->node_type==PAIR_FRONT_NODE);

        if((it->node_type)==PAIR_END_NODE)
        {
            //shc_log_info(shc_logname, "p1 end, node id %u -> node id %u\n",
            //        graph[vd_old].node_id, graph[vd_new].node_id);

            p1_end[it->align_read_id] = vd_new;
            graph[vd_new].term_nodes_info.emplace_back(PAIR_END_NODE, it->align_read_id);
        }
        else if((it->node_type)==PAIR_FRONT_NODE)
        {
            //shc_log_info(shc_logname, "p2 front, node id %u -> node id %u\n",
            //        graph[vd_old].node_id, graph[vd_new].node_id);

            p2_front[it->align_read_id] = vd_new;
            graph[vd_new].term_nodes_info.emplace_back(PAIR_FRONT_NODE, it->align_read_id);
        }
        else
        {
            shc_log_error("cannot transfer NO PAIR TYPE READ\n");
            exit(1);
        }
    }

    v1_term_nodes_info.clear();
}

void Sequence_graph_handler::condense_graph()
{
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start condense graph\n");
#endif

    start_timer(&timer);
    out_ei_t out_ei, out_ei_end;
    vi_t vi, vi_end, next;
    boost::tie(vi, vi_end) = boost::vertices(graph);
    //walk through the graph list, and update if
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "before condense, number of vertices is %d, number of edge is %d\n",
            boost::num_vertices(graph), boost::num_edges(graph));
#endif
    std::set<vd_t> vd_remove_set;
    uint64_t num_vertices_removed = 0;
    for(next=vi; vi!=vi_end; vi=next)
    {
        //shc_log_info(shc_logname, "next node\n");
        ++next;
        vd_t vd_source = *vi;
        std::set<vd_t> local_rm_set;
        //shc_log_info(shc_logname, "Ori: %s \n", graph[vd_source].seq.c_str());


        vd_t vd_new = local_condense(vd_source, local_rm_set);
        if(vd_new != NULL) // something condensed
        {
            //shc_log_info(shc_logname, "*************new node************* \n");
            //log_node_info(vd_new, true, false, false);
            //shc_log_info(shc_logname, "condensed %s \n", graph[vd_new].seq.c_str());
            for(std::set<vd_t>::iterator it=local_rm_set.begin();
                it!=local_rm_set.end(); ++it)
            {
                assert(boost::in_degree(*it, graph)==0 &&
                        boost::out_degree(*it, graph)==0);
                //shc_log_info(shc_logname, "************ condensed nodes****\n");
                //log_node_info(*it, false, false);

                vd_remove_set.insert(*it);
            }
        }
        //else
        //{
        //    shc_log_info(shc_logname, "nothing\n");
        //}
    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "before remove nodes\n");
#endif
    for(std::set<vd_t>::iterator vrs_it=vd_remove_set.begin();
                    vrs_it != vd_remove_set.end(); vrs_it++)
    {
        //shc_log_info(shc_logname, "in %d, out %d\n", boost::in_degree(*it, graph),
        //log_node_info(*vrs_it, false, false);
        //shc_log_info(shc_logname, "before remove nodes\n");
        assert(1==exist_vd.erase(*vrs_it));
        boost::remove_vertex(*vrs_it, graph);
        //shc_log_info(shc_logname, "after remove nodes\n");
    }
    //find_all_simple_pair();

    //shc_log_info(shc_logname, "Finish condense graph\n");
    //for(std::set<vd_t>::iterator it=vd_remove_set.begin(); it != vd_remove_set.end(); it++)
    //    boost::remove_vertex(*it,graph);
#ifdef LOG_SEQ_GRAPH
shc_log_info(shc_logname, "after condense, number of vertices is %d, number of edge is %d\n",
        boost::num_vertices(graph), boost::num_edges(graph));
    shc_log_info(shc_logname, "number vertices removed is %u\n",
                                                        num_vertices_removed);

#endif
#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "condense graph and remove nodes finish ";
    stop_timer(&timer);
#endif
}

int Sequence_graph_handler::get_num_xnodes()
{
    int num_xnode = 0;
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        if(is_xnode(*it))
        {
            num_xnode++;
        }
    }
    return num_xnode;
}

void Sequence_graph_handler::update_xnode_set()
{
    start_timer(&timer);
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
    //shc_log_info(shc_logname, "Finish getting all xnodes\n");
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
    Read_acc acc = coll_read_list.get_read(read_info.read_id);
    char left_read_base = *(acc.read_ptr + read_info.start - 1);
    in_eip_t in_eip = boost::in_edges(vd, graph);
    for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
    {
        vd_t vd_source = boost::source(*in_ei, graph);
        bundled_node_p & source_node = graph[vd_source];
        bundled_edge_p & in_edge = graph[*in_ei];

        if(source_node.seq_len()==in_edge.weight)
        {
#ifdef LOG_SEQ_GRAPH
            shc_log_info(shc_logname ,"encounter same length,check read bridge left\n");
#endif
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
    Read_acc acc = coll_read_list.get_read(read_info.read_id);
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
#ifdef LOG_SEQ_GRAPH
            shc_log_info(shc_logname ,"encounter same length,check read bridge right\n");
#endif
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
    std::vector<Bdg_read_info> & reads_info = graph[vd].reads_info;

    std::vector<Bdg_read_info> valid_reads_info;

    for(std::vector<Bdg_read_info>::iterator it = reads_info.begin();
                                             it!= reads_info.end(); it++)
    {
        //shc_log_info(shc_logname, "read id %u\n", it->read_id);
        Read_acc acc = coll_read_list.get_read(it->read_id);

        if(!is_read_bridge_node(acc.read_ptr, acc.len, it->start, vd))
            continue;

        if(is_read_bridge_left(*it, vd) && is_read_bridge_right(*it,vd))
        {
            valid_reads_info.push_back(*it);
        }
    }

    reads_info.assign(valid_reads_info.begin(), valid_reads_info.end());
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
        Read_acc acc = coll_read_list.get_read(ri_it->read_id);

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
        if(in_base_set.size()==4 && out_base_set.size()==4)
            break;
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

uint64_t Sequence_graph_handler::find_all_simple_pair()
{
    //log_all_node(true, false);
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start find_all_simple_pair\n");
#endif
    uint64_t num_pair = p1_end.size();
    assert(p1_end.size() == p2_front.size());
    uint64_t num_useful_sr = 0;
    int64_t num_sr = 0;
    for(uint64_t i=0; i<num_pair; i++)
    {
        bool is_sr =false;
        if(check_and_add_simple_read_pair(i, is_sr))
            num_useful_sr++;
        if(is_sr)
            num_sr++;

    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "find %u useful super read\n", num_useful_sr);
    shc_log_info(shc_logname, "Finish find_all_simple_pair\n");
#endif
    return num_useful_sr;
}


void Sequence_graph_handler::clear_term_info(uint64_t align_index)
{
    vd_t vd1 = p1_end[align_index];
    vd_t vd2 = p2_front[align_index];
    graph[vd1].term_nodes_info.clear();
    graph[vd2].term_nodes_info.clear();
    p1_end[align_index] = NULL;
    p2_front[align_index] = NULL;
}

/**
 * This function checks if the pair read is simply connected.
 * '_' represents reads. The following chart illustrate a simple connect read
 * ____1          2_____        *                 __2__
 * xxxx            xxxxx        * xxxx            xxxxx
 *      \ _1_  _2 /             *      \  _1_  _2_ /
 *        xxxxxxx               *        xxxxxxxxx
 *      /         \             *      /           \
 * xxxxx           xxxxxx       * xxxxx             xxxxxx
 *  useful for multibridge      * just recognize full pair
 * Input : alignment read index
 * return : True if super read is formed and useful for multi-bridge
 *          False if no super read is formed or if the super read is not useful
 *          for multibridge purpose
 */
bool Sequence_graph_handler::check_and_add_simple_read_pair(uint64_t i, bool & is_sr)
{
    //shc_log_info(shc_logname, "start\n");
    bool is_multibridge_useful = false;
    vd_t vd1 = p1_end[i];
    vd_t vd2 = p2_front[i];
    if(vd1==NULL || vd2==NULL)
    {
        //shc_log_info(shc_logname, "one of node is NULL\n");
        clear_term_info(i);
        return false;
    }
    if( vd1 == vd2)
    {
#ifdef LOG_PAIR_SEARCH
        shc_log_info(shc_logname, "seq v1 node %u\n", graph[vd1].node_id);
        shc_log_info(shc_logname, "seq v2 node %u\n", graph[vd2].node_id);
        shc_log_info(shc_logname, "node: %s\n", graph[vd2].seq.c_str());
#endif
        std::string middle;

        read_num_t l_start_read_id = coll_read_list.convert_align_index_to_local_index(i, true);
        read_num_t l_stop_read_id = coll_read_list.convert_align_index_to_local_index(i, false);

        Read_acc p1_acc = coll_read_list.p1_reads.get_read(l_start_read_id);
        log_read_seq(p1_acc);
        Read_acc p2_acc = coll_read_list.p2_reads.get_read(l_stop_read_id);
        log_read_seq(p2_acc);

        read_num_t v_start_read_id =
                coll_read_list.convert_local_index_to_virtual_index(l_start_read_id, true);
        read_num_t v_stop_read_id =
                coll_read_list.convert_local_index_to_virtual_index(l_stop_read_id, false);

        curr_align_index = i;
        int64_t middle_len;
        bool is_mismatch = false;;
        if(find_middle_seq(vd1, l_start_read_id, l_stop_read_id, middle, middle_len, is_mismatch))
        {
            is_multibridge_useful = true;
#ifdef LOG_PAIR_SEARCH
            shc_log_info(shc_logname, "useful\n");
            shc_log_warning("useful\n");
#endif
        }


        if(is_mismatch)
        {
            shc_log_info(shc_logname, "find a mismatch\n");
            clear_term_info(i);
            is_sr = false;
            return false;
        }
#ifdef LOG_PAIR_SEARCH
        shc_log_info(shc_logname, "middle: %s\n",middle.c_str());
#endif
        is_sr = true;
        if(middle_len<=0)
        {
            shc_log_info(shc_logname, "find an overlap pair read\n");
        }
        coll_read_list.add_super_read(v_start_read_id, v_stop_read_id,
                                      l_start_read_id, l_stop_read_id, middle);

        //debug
#ifdef LOG_PAIR_SEARCH
        Super_read sr(l_start_read_id, l_stop_read_id, middle);
        Read_acc acc = coll_read_list.get_super_read(sr);
        std::string sr_str(acc.read_ptr, acc.len);
        shc_log_info(shc_logname, "find simple sr %s\n", sr_str.c_str());
#endif
        clear_term_info(i);
        if(is_multibridge_useful)
        {
            return true;
        }
    }
    //shc_log_info(shc_logname, "end\n");
    return false;
}

int64_t Sequence_graph_handler::
find_terminal_index(vd_t vd, std::string & seq, char * read_ptr, read_length_t len,
                                             read_length_t acc_len, bool is_p1)
{
    if(len >= acc_len)
        return -1;

    if(is_p1)
    {
        std::string read_p1_last(read_ptr, len);
        int64_t p1_end_index = seq.rfind(read_p1_last)+kmer_length; //  where p1 last kmer kmer located on node seq
        int64_t p1_end_index_h = seq.find(read_p1_last)+kmer_length;

        if(p1_end_index==std::string::npos )
        {
            shc_log_info(shc_logname, "is super read %s\n",
                    (coll_read_list.is_sr(curr_align_index))?("yes"):("no") );
            shc_log_info(shc_logname, "Index: %u, Node ID %u: seq %s cannot find p1 last kemr %s\n",
                curr_align_index, graph[vd].node_id, seq.c_str(), read_p1_last.c_str());
            shc_log_error("comp %d, do not find p1 last kmer in the node\n", curr_comp);
            exit(1);
        }

        if(p1_end_index!=p1_end_index_h)
        {
            return find_terminal_index(vd, seq, read_ptr-1, len+1, acc_len, is_p1);
        }
        else
        {
            return p1_end_index;
        }
    }
    else
    {
        std::string read_p2_first(read_ptr, len);
        int64_t p2_first_index = seq.find(read_p2_first);
        int64_t p2_first_index_h = seq.rfind(read_p2_first);
        if(p2_first_index==std::string::npos)
        {
            shc_log_info(shc_logname, "is super read %s\n",
                    (coll_read_list.is_sr(curr_align_index))?("yes"):("no") );
            shc_log_info(shc_logname, "Index: %u, Node ID %u: seq %s cannot find p2 front kemr %s\n",
                curr_align_index, graph[vd].node_id, seq.c_str(), read_p2_first.c_str());
            shc_log_error("comp %d, do not find p2 first kmer in the node\n", curr_comp);
            exit(1);
        }

        if(p2_first_index!=p2_first_index_h)
        {
            return find_terminal_index(vd,seq, read_ptr, len+1, acc_len, is_p1);
        }
        else
        {
            return p2_first_index;
        }
    }


}

bool Sequence_graph_handler::find_middle_seq(vd_t vd, read_num_t local_read_1,
                read_num_t local_read_2, std::string & middle, int64_t & middle_seq_len,
                bool & is_mismatch)
{
    bundled_node_p & curr_node = graph[vd];
    // get node seq
    std::string & seq = curr_node.seq;

    // get read terminal sub seq of kmer length
    Read_acc p1_acc = coll_read_list.p1_reads.get_read(local_read_1);
    Read_acc p2_acc = coll_read_list.p2_reads.get_read(local_read_2);

    // get index
    int64_t p1_end_index  = find_terminal_index( vd,
            seq, p1_acc.read_ptr+p1_acc.len-kmer_length, kmer_length, p1_acc.len, true);
    int64_t p1_first_index = p1_end_index - p1_acc.len;

    int64_t p2_first_index =  find_terminal_index(vd,
                             seq, p2_acc.read_ptr, kmer_length, p2_acc.len, false);
    int64_t p2_end_index = p2_first_index + p2_acc.len;

    if(p1_end_index < 0 || p2_first_index < 0 ||
       p1_end_index >= p2_end_index || p2_first_index <=  p1_first_index)
    {
        shc_log_info(shc_logname, "find mismatch\n");
        is_mismatch = true;
        return false;
    }


    //shc_log_info(shc_logname, "p1_last_index %u\n", p1_end_index);
    //shc_log_info(shc_logname, "p2_first_index %u\n", p2_first_index);

    middle_seq_len = p2_first_index - p1_end_index;
    if(middle_seq_len > 0)
    {
        middle = seq.substr(p1_end_index, middle_seq_len);
    }
    else if(middle_seq_len==0)
    {
        middle = "";
    }
    else
    {
        middle = OVERLAP_CHAR + std::to_string(-middle_seq_len);
    }

    //determine if super read useful
    if ( p1_end_index>(p1_acc.len-kmer_length)||
        (seq.size()-p2_first_index) > p2_acc.len )
    {
        shc_log_info(shc_logname, "seq fully cover one of the read or both\n");
        return false;   // seq fully cover one of the read
    }
    else
    {
        //shc_log_info(shc_logname, "end useful\n");
        return true;
    }
}

void Sequence_graph_handler::
get_xnodes(std::set<vd_t> & new_nodes, std::set<vd_t> & xnodes)
{
    //shc_log_info(shc_logname, "Start %d new nodes\n", new_nodes.size());
    if(xnodes.size()!=0)
        xnodes.clear();
    //walk through the graph list, and update if
    for(std::set<vd_t>::iterator it=new_nodes.begin(); it!=new_nodes.end(); it++)
    {
        if(is_xnode(*it))
        {
            xnodes.insert(*it);
        }
    }
    //shc_log_info(shc_logname, "End %d xnodes\n", xnodes.size());
}


void Sequence_graph_handler::bridge_all_xnodes()
{

    start_timer(&timer);

#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start bridging xnode\n");
#endif
    int num_iter = 0;

    std::set<vd_t> nodes_visited;
    //shc_log_info(shc_logname, "\t\tbefore resolve\n");
    //log_term_array(false);

    start_timer(&resolve_pair_timer);
    if(setting.has_pair && get_num_edges()>0 &&
        setting.seq_graph_setup.max_hop_pair >= 0)
    {
        std::set<vd_t> vd_to_be_rechecked;
        if(resolve_all_pair_reads(setting.seq_graph_setup.max_hop_pair,
          vd_to_be_rechecked) > 0)
        {
            for(std::set<vd_t>::iterator it=vd_to_be_rechecked.begin();
                                        it!=vd_to_be_rechecked.end(); it++)
            {
                std::set<vd_t>::iterator n_it = nodes_visited.find(*it);
                if(n_it != nodes_visited.end())
                    nodes_visited.erase(n_it);
            }
        }
    }
    #ifdef PRINT_TIME
        std::cout << " finish resolve_pair_timer" << std::endl;
        stop_timer(&resolve_pair_timer);
    #endif
    //std::cout << "Finish resolve " << std::endl;
    //stop_timer(&resolve_pair_time);
    //shc_log_info(shc_logname, "finsih resolve pair\n");
    stop_timer_np(&resolve_pair_timer);
    //shc_log_info(shc_logname, "\t\tafter resolve\n");
    //log_term_array(false);

    //std::ofstream f_write(setting.local_files.output_path + "/analysis/C_after_condense");
    mb_condense_timer.time_us = 0;
    while(true)
    {
#ifdef PRINT_TIME
        start_timer(&test_timer);
        int num_read = 0;
#endif
        num_edge_added = 0;
        int num_xnode_bridged = 0;
        //shc_log_info(shc_logname, "start new iter %d\n", num_iter);

        //std::cout << "num nodes "<< get_num_nodes() << std::endl;
        //std::cout << "num edges "<< get_num_edges() << std::endl;
        //std::cout << "num isolated nodes " << get_num_isolated_nodes() << std::endl;
        //std::cout << "num xnode " << xnode_set.size() << std::endl;
        //if (num_iter >= 2)
        //    exit(0);
        //log_graph_to_file(graph_writer, num_iter);

        for(Node_set_iterator x_iter=xnode_set.begin();
                            x_iter!=xnode_set.end(); ++x_iter)
        {
            //if(nodes_visited.find(*x_iter) == nodes_visited.end())
            //{
                nodes_visited.insert(*x_iter);
                bundled_node_p & node = graph[*x_iter];
                refresh_bridging_reads(*x_iter);
                if(is_xnode_bridged(*x_iter))
                {
#ifdef PRINT_TIME
                    num_read += graph[*x_iter].reads_info.size();
#endif
                    perform_bridging(*x_iter);
                    num_xnode_bridged++;
                }
            //}
        }
        total_num_bridged += num_xnode_bridged;

#ifdef LOG_SEQ_GRAPH
        std::cout << "comp " << curr_comp << " num iter " << num_iter
                  << " num_xnode_bridged " << num_xnode_bridged
                  << " nodes_visited " << nodes_visited.size()  << std::endl;
        shc_log_info(shc_logname, "comp %d, num iter %d num_xnode_bridged %d\n",
                        curr_comp, num_iter, num_xnode_bridged);
#endif



#ifdef PRINT_TIME
        std::cout << "comp " << curr_comp << " num iter " << num_iter
                  << " num_xnode_bridged " << num_xnode_bridged
                  << " nodes_visited " << nodes_visited.size() << ". Checked "
                  << num_read << " reads"
                  << "mb_condense_timer takes " << mb_condense_timer.time_us << " us"
                  << "num node " << get_num_nodes() << ", num edge " << get_num_edges()
                  << ", num_edge_added " << num_edge_added << std::endl;

        mb_condense_timer.time_us = 0;
        shc_log_info(shc_logname, "comp %d, num iter %d num_xnode_bridged %d\n",
                        curr_comp, num_iter, num_xnode_bridged);
        stop_timer(&test_timer);
#endif

        //condense_graph();

        num_iter++;
        if(num_xnode_bridged==0)
        {
            condense_graph();
            return;
        }
        else
        {
            update_xnode_set();
        }
    }

}

void Sequence_graph_handler::perform_bridging(vd_t vd)
{
    bundled_node_p & curr_node = graph[vd];
    //shc_log_info(shc_logname, "Node %u\n", curr_node.node_id);
    assert(curr_node.reads_info.size() > 0);
    assert(boost::in_degree(vd, graph)>=2 && boost::out_degree(vd, graph)>=2);

    //shc_log_info(shc_logname, "\n");
    std::vector<vd_t> u_list, w_list;
    //std::cout << "bridging " << curr_node.seq << std::endl;
    //shc_log_info(shc_logname, "bridging %s\n", curr_node.seq.c_str());

    bool contain_loop = false;
    vd_t u_loop_node, w_loop_node;
    edge_weight_t loop_weight;
    double loop_count;

    //log_node_info(vd, true, false, false);

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
            loop_count = graph[*it].count;
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
    shc_log_info(shc_logname, "left extension %d\n", u_list.size());
    for(std::vector<vd_t>::iterator it=u_list.begin(); it!=u_list.end(); it++)
    {
        shc_log_info(shc_logname, "%s\n", graph[*it].seq.c_str());
    }

    shc_log_info(shc_logname, "right extension %d\n", w_list.size());
    for(std::vector<vd_t>::iterator it=w_list.begin(); it!=w_list.end(); it++)
    {
        shc_log_info(shc_logname, "%s\n", graph[*it].seq.c_str());
    }
#endif


    // deal with self loop
    if (contain_loop)
    {
#ifdef LOG_SEQ_GRAPH
        shc_log_info(shc_logname, "contain self loop\n");
#endif
        aer_t aer = boost::add_edge(w_loop_node, u_loop_node, graph);
        graph[aer.first] = bundled_edge_p(loop_weight+2, loop_count);
        num_edge_added++;
        //shc_log_info(shc_logname, "self loop\n");
        //shc_log_info(shc_logname, "%s link to %s\n", graph[w_loop_node].seq.c_str(), graph[u_loop_node].seq.c_str());
    }

    // clear edges around v, and vertex will be destroyed at the end
    boost::clear_vertex(vd, graph);

    // link extension nodes
    std::vector<Bdg_read_info> & reads_info = curr_node.reads_info;

    std::vector<vd_t> matched_u;
    std::vector<vd_t> matched_w;

    std::set<vd_t> bridged_u;
    std::set<vd_t> bridged_w;

    // make sure reads of vertex are not repeated
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "bridging between extension nodes\n");
#endif
    //shc_log_info(shc_logname, "start read\n");
    for(std::vector<Bdg_read_info>::iterator ri_it=reads_info.begin();
                                      ri_it!=reads_info.end(); ++ri_it)
    {
        // if the same read is processed, ignore the rest
        Read_acc acc = coll_read_list.get_read(ri_it->read_id);
        //std::string read(acc.read_ptr, acc.len);
        //shc_log_info(shc_logname, "read unique id %d, %s\n", unique_read_index, read.c_str());
        //shc_log_info(shc_logname, "bridge u size %d, bridge w size %d\n", bridged_u.size(), bridged_w.size());

        char * l_seq_start = acc.read_ptr + ri_it->start - 1;
        for(std::vector<vd_t>::iterator
                node_it=u_list.begin(); node_it!=u_list.end(); ++node_it)
        {
            bundled_node_p & node_p = graph[*node_it];
            //shc_log_info(shc_logname, "left ext %s\n", node_p.seq.c_str());
            if(node_p.seq.at(0) == *l_seq_start)//memcmp(&node_p.seq.at(0), l_seq_start, node_p.seq_len())==0)
            {
               // shc_log_info(shc_logname, "matched left ext %s\n",node_p.seq.c_str());
                matched_u.push_back(*node_it);
                bridged_u.insert(*node_it);
            }
        }

        char * seq_start = acc.read_ptr + ri_it->start;
        for(std::vector<vd_t>::iterator
                node_it=w_list.begin(); node_it!=w_list.end(); ++node_it)
        {
            bundled_node_p & node_p = graph[*node_it];
            //shc_log_info(shc_logname, "right ext %s\n", node_p.seq.c_str());
            if(node_p.seq.back()== *(seq_start+node_p.seq.size()-1))//memcmp(&node_p.seq.at(0), seq_start, node_p.seq_len())==0)
            {
                //shc_log_info(shc_logname, "matched right ext %s\n", node_p.seq.c_str());
                matched_w.push_back(*node_it);
                bridged_w.insert(*node_it);
            }
        }
        if(matched_u.size()!=1 || matched_w.size()!=1)
        {

            shc_log_info(shc_logname, "in special case that some kmer in read is missing while building the graph\n");

            matched_u.clear();
            matched_w.clear();
            continue;
        }

        // create edge between nodes and link read
        vd_t u = matched_u[0];
        vd_t w = matched_w[0];
        node_link_read(u, ri_it->read_id, ri_it->start-1);
        node_link_read(w, ri_it->read_id, ri_it->start);

        //if(graph[vd].seq == "CCCCCAGCCGCAGGAGCCGAACCCCCAGCCG")
        //{
        //    shc_log_info(shc_logname, "special link from %c to %c\n", graph[u].seq.front(), graph[w].seq.back());
        //    log_read_seq(coll_read_list.get_read(ri_it->read_id));
        //}

        //if there is no such edge
        aer_t aer = boost::edge(u,w,graph);
        if(!aer.second)
        {
            //shc_log_info(shc_logname, "%s link to %s\n", graph[u].seq.c_str(), graph[w].seq.c_str());
            aer = boost::add_edge(u, w, graph);
            graph[aer.first] = bundled_edge_p(curr_node.seq_len(), acc.read_count);
            num_edge_added++;
        }
        else
        {
            graph[aer.first].count+=acc.read_count;
        }
        matched_u.clear();
        matched_w.clear();
    }
    //shc_log_info(shc_logname, "finish read\n");

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
        graph[aer.first] = bundled_edge_p(curr_node.seq_len(), 1);
        num_edge_added++;

        //shc_log_info(shc_logname, "special %s link to %s\n", graph[*u_it].seq.c_str(), graph[*w_it].seq.c_str());

#ifdef LOG_SEQ_GRAPH
        std::cout << "speical bridge" << std::endl;
        shc_log_info(shc_logname, "special bridging apply\n");
#endif
    }
    else
    {
        assert(u_list.size()-bridged_u.size()==0 &&
               w_list.size()-bridged_w.size()==0);

#ifdef LOG_SEQ_GRAPH
        std::cout << "full bridge" << std::endl;
        shc_log_info(shc_logname, "everything directly bridged\n");
#endif
    }

    //condense graph
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "local condense bridged nodes\n");
#endif

    //shc_log_info(shc_logname, "\t\tAfter bridge\n");

    local_condense_bridged(u_list, w_list);


    boost::remove_vertex(vd, graph);
    assert(1==exist_vd.erase(vd));
    /*
    for(int new_vd_i = 0; new_vd_i < w_list.size(); new_vd_i++)
    {
        shc_log_info(shc_logname, "new node %d\n", new_vd_i);
        log_node_info(w_list[new_vd_i], true, true, false);
    }
    shc_log_info(shc_logname, "\t\tAfter bridge end\n");
    */
}

// This function takes one argument, condense left or right if it can,
// otherwise return false
bool Sequence_graph_handler::
condense(vd_t vd, bool is_check_left , vd_t vd_new, Read_set & read_set,
                                        int & seq_len_before_last_concat)
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
            new_node.count += source_node.count;
            new_node.prevalence += source_node.prevalence;
            bundled_edge_p & edge_p = graph[*(in_eip.first)];
            //shc_log_info(shc_logname, "condense left\n");
            transfer_read_align_index(vd_source, vd_new);

            std::string seq_back;
            seq_back.resize(source_node.seq.size());
            std::reverse_copy(source_node.seq.begin(), source_node.seq.end(), seq_back.begin());
            seq_len_before_last_concat = new_node.seq_len();
            new_node.seq += seq_back.substr(edge_p.weight);

            //update previous reads info
            /*
            std::vector<Bdg_read_info> new_set;
            for(Read_set_iterator it = read_set.begin(); it!=read_set.end(); it++)
            {
                if (it->start > source_node.seq_len() - edge_p.weight)
                {
                    read_length_t new_start = it->start-source_node.seq_len()+edge_p.weight;
                    Bdg_read_info read_info(it->read_id, new_start);
                    new_set.push_back(read_info);
                }
            }
            // add new reads
            read_set.clear();
            read_set.insert(new_set.begin(), new_set.end());
            read_set.insert(source_node.reads_info.begin(), source_node.reads_info.end());
            */
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
            new_node.count += target_node.count;
            new_node.prevalence += target_node.prevalence;
            // info is only updated in the left side
            bundled_edge_p & edge_p = graph[*(out_eip.first)];

            //shc_log_info(shc_logname, "condense right\n");
            transfer_read_align_index(vd_target, vd_new);

            seq_len_before_last_concat = new_node.seq_len();
            new_node.seq = new_node.seq + target_node.seq.substr(edge_p.weight);

            /*
            //read_set.insert(target_node.reads_info.begin(), target_node.reads_info.end());
            for(std::vector<Bdg_read_info>::iterator it = target_node.reads_info.begin();
                        it!=target_node.reads_info.end(); it++)
            {
                if (it->start > seq_len_before_concat - edge_p.weight)
                {
                    read_length_t new_start = it->start-seq_len_before_concat+edge_p.weight;
                    Bdg_read_info read_info(it->read_id, new_start);
                    read_set.insert(read_info);
                }
            }
            */
        }
    }

    assert(boost::out_degree(vd_source, graph) == 1 &&
           boost::in_degree(vd_target, graph) == 1);

    return true;
}

void Sequence_graph_handler::special_condense_self_loop()
{
    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            if(vd == vd_target)
            {
                bundled_node_p & curr_node = graph[vd];
                bundled_edge_p & edge_p = graph[*out_ei];
                std::string new_seq = curr_node.seq + curr_node.seq.substr(edge_p.weight);
                curr_node.seq = new_seq;
            }
        }
    }
}

void Sequence_graph_handler::
condense_left(vd_t vd_target, vd_t vd_new, std::set<vd_t> & vd_remove_set,
                                           Read_set & read_set)
{
    vd_t vd_source;
    int seq_len_before_last_concat = 0;
    while (condense(vd_target, true, vd_new, read_set, seq_len_before_last_concat))
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
    read_set.insert(graph[vd_leftmost].reads_info.begin(), graph[vd_leftmost].reads_info.end());

    //link left nodes to new nodes
    in_eip_t in_eip = boost::in_edges(vd_leftmost, graph);
    for(in_ei_t it = in_eip.first; it!=in_eip.second; it++)
    {
        aer_t aer = boost::add_edge(boost::source(*it, graph), vd_new, graph);
        //copy the edge count
        graph[aer.first] = bundled_edge_p(graph[*it].weight, graph[*it].count);
    }
    boost::clear_in_edges(vd_leftmost, graph);
}


void Sequence_graph_handler::
condense_right(vd_t vd_source, vd_t vd_new ,std::set<vd_t> & vd_remove_set,
                                           Read_set & read_set )
{
    vd_t vd_target;
    int seq_len_before_last_concat = 0;
    int last_edge_weight = 0;
    while (condense(vd_source, false, vd_new, read_set, seq_len_before_last_concat))
    {
        //get the source node and clear
        out_eip_t out_eip = boost::out_edges(vd_source, graph);
        vd_target = boost::target(*(out_eip.first), graph);
        last_edge_weight = graph[*(out_eip.first)].weight;
        boost::clear_out_edges(vd_source, graph);
        vd_remove_set.insert(vd_target);


        //prepare for next iteration
        vd_source = vd_target;
    }

    vd_t vd_rightmost = vd_source; // the node which stands for the right end

    std::vector<Bdg_read_info> & reads_info = graph[vd_rightmost].reads_info;
    for(std::vector<Bdg_read_info>::iterator it=reads_info.begin();
                                it!=reads_info.end(); it++)
    {

        if (static_cast<int>(it->start) > seq_len_before_last_concat -
                                    last_edge_weight)
        {
            read_length_t new_start = it->start - seq_len_before_last_concat
                                            + last_edge_weight;
            Bdg_read_info read_info(it->read_id, new_start);

            read_set.insert(read_info);
        }

    }

    //link right node
    out_eip_t out_eip = boost::out_edges(vd_rightmost, graph);
    for(out_ei_t it = out_eip.first; it!=out_eip.second; it++)
    {
        aer_t aer = boost::add_edge(vd_new, boost::target(*it, graph), graph);
        //copy the edge
        graph[aer.first] = bundled_edge_p(graph[*it].weight, graph[*it].count);
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

bool Sequence_graph_handler::is_only_self_loop(vd_t vd)
{
    if(boost::in_degree(vd, graph)!= 1)
        return false;
    if(boost::out_degree(vd, graph)!= 1)
        return false;
    out_eip_t out_eip = boost::out_edges(vd, graph);
    vd_t vd_target = boost::target(*(out_eip.first), graph);
    if(vd_target == vd)
        return true;
}

bool Sequence_graph_handler::is_contain_self_loop(vd_t vd, ed_t & self_edge)
{
    out_eip_t out_eip = boost::out_edges(vd,graph);
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        if(vd_target == vd)
        {
            self_edge = *out_ei;
            return true;
        }
    }
    return false;
}

//return ture when special case applies
bool Sequence_graph_handler::py_special_condense(vd_t vd)
{
    ed_t self_edge;
    if(is_contain_self_loop(vd, self_edge))
    {
        bundled_node_p & curr_node = graph[vd];
        bundled_edge_p & curr_edge = graph[self_edge];
        std::string new_seq = curr_node.seq + curr_node.seq.substr(curr_edge.weight);
        graph[vd].seq = new_seq;
        boost::remove_edge(self_edge, graph);
        return true;
    }
    else
    {
        return false;
    }
}

// always remember to clear the vd_remove_set afterward
vd_t Sequence_graph_handler::
local_condense(vd_t vd, std::set<vd_t> & vd_remove_set)
{
    //py_special_condense(vd);

    bool left_condensable = is_condensable(vd, true);
    bool right_condensable = is_condensable(vd, false);

    Read_set read_coll; // in case of repeat

    read_coll.set_empty_key(read_empty);
    read_coll.set_deleted_key(read_delete);
    if( left_condensable || right_condensable)
    {
        vd_t vd_new = boost::add_vertex(graph);
        exist_vd.insert(vd_new);
        graph[vd_new] = bundled_node_p(graph[vd].seq, glob_node_id++);

        graph[vd_new].reads_info.clear();
        graph[vd_new].term_nodes_info.clear();
        transfer_read_align_index(vd, vd_new);
        //accurate_start_timer(&mb_condense_timer);
        //read_coll.insert(graph[vd].reads_info.begin(), graph[vd].reads_info.end());
        //accurate_accumulate_timer(&mb_condense_timer);
        //log_node_info(*n_it, true, false);
        //shc_log_info(shc_logname, "before condense left\n");
        std::reverse(graph[vd_new].seq.begin(), graph[vd_new].seq.end());
        condense_left(vd, vd_new, vd_remove_set, read_coll);

        std::reverse(graph[vd_new].seq.begin(), graph[vd_new].seq.end());
        //shc_log_info(shc_logname, "before condense right\n");
        condense_right(vd, vd_new ,vd_remove_set, read_coll);
        //std::cout << "read_coll size " << read_coll.size() << std::endl;
        for(Read_set_iterator it= read_coll.begin();
                                              it!=read_coll.end(); it++)
        {
            graph[vd_new].reads_info.push_back(*it);
        }
        //graph[vd_new].reads_info = graph[vd].reads_info;
        //refresh_bridging_reads(vd_new);
        boost::clear_vertex(vd, graph);
        vd_remove_set.insert(vd);

        /*
            std::cout << "before" << std::endl;
            std::cout << graph[vd_new].seq << std::endl;
            shc_log_info(shc_logname, "seq %s\n", graph[vd_new].seq.c_str());
            for(std::vector<Terminal_node_info>::iterator
                    it= graph[vd_new].term_nodes_info.begin();
                    it!=graph[vd_new].term_nodes_info.end(); it++)
            {
                shc_log_info(shc_logname, "align_I %u\n", it->align_read_id);

                Read_acc acc =  coll_read_list.get_read(
                    coll_read_list.convert_align_index_to_virtual_index(it->align_read_id, it->node_type));
                std::string read(acc.read_ptr, acc.len);
                shc_log_info(shc_logname, "%s, index %u, %s\n",
                            ((it->node_type==PAIR_END_NODE)?("PAIR_END_NODE"):("PAIR_FRONT_NODE")),
                            it->align_read_id, read.c_str());
                std::cout << ((it->node_type==PAIR_END_NODE)?("PAIR_END_NODE"):("PAIR_FRONT_NODE"))
                          << "  index " << it->align_read_id  << std::endl;

            }
            is_node_contain_terminals(vd_new);
        }
        */

        //is_node_contain_terminals(vd_new);
        /*
        if(graph[vd_new].node_id==270784)
        {
            std::cout << "after" << std::endl;
            std::cout << graph[vd_new].seq << std::endl;
            for(std::vector<Terminal_node_info>::iterator
                    it= graph[vd_new].term_nodes_info.begin();
                    it!=graph[vd_new].term_nodes_info.end(); it++)
            {
                std::cout << ((it->node_type==PAIR_END_NODE)?("PAIR_END_NODE"):("PAIR_FRONT_NODE"))
                          << "  index " << it->align_read_id  << std::endl;
            }
        }
        */

        return vd_new;
    }
    else
    {
        return NULL;//shc_log_info(shc_logname, "condense nothing\n");
    }
}

bool Sequence_graph_handler::is_node_contain_terminals(vd_t vd)
{
    std::vector<Terminal_node_info> & term_nodes_info =
                                            graph[vd].term_nodes_info;
    std::string & seq = graph[vd].seq;
    //assert(!v1_term_nodes_info.empty());
    for(std::vector<Terminal_node_info>::iterator
        it=term_nodes_info.begin(); it != term_nodes_info.end(); it++)
    {

        if(it->node_type == PAIR_END_NODE)
        {
            read_num_t local_read_1 =
                    coll_read_list.convert_align_index_to_local_index(it->align_read_id, true);
            Read_acc p1_acc = coll_read_list.p1_reads.get_read(local_read_1);
            int64_t p1_end_index  = find_terminal_index(vd,
                seq, p1_acc.read_ptr+p1_acc.len-kmer_length, kmer_length, p1_acc.len, true);

            //shc_log_info(shc_logname, "p1 last %s\n", seq.substr(p1_end_index-kmer_length).c_str());
            if(p1_end_index < 0 )
            {
                shc_log_error("comp %d cannot find\n",  curr_comp);
                exit(1);
            }
        }
        else if (it->node_type == PAIR_FRONT_NODE)
        {
            read_num_t local_read_2 =
                    coll_read_list.convert_align_index_to_local_index(it->align_read_id, false);
            Read_acc p2_acc = coll_read_list.p2_reads.get_read(local_read_2);
            int64_t p2_first_index =  find_terminal_index(vd,
                                 seq, p2_acc.read_ptr, kmer_length, p2_acc.len, false);
            //shc_log_info(shc_logname, "p2 front %s\n", seq.substr(p2_first_index, kmer_length).c_str());
            if(p2_first_index <0)
            {
                shc_log_error("comp %d cannot find\n",  curr_comp);
                exit(1);
            }
        }
        else
        {
            shc_log_error("comp %d no type\n",  curr_comp);
            exit(1);
        }
    }
    return true;
}



// all nodes have been bridged
void Sequence_graph_handler::
local_condense_bridged(std::vector<vd_t> & u_list, std::vector<vd_t> & w_list)
{
    //shc_log_info(shc_logname, "Start local condense\n");
    std::set<vd_t> vd_remove_set;
    u_list.insert(u_list.end(), w_list.begin(), w_list.end());
    w_list.clear();

    for(std::vector<vd_t>::iterator n_it=u_list.begin(); n_it!=u_list.end(); ++n_it)
    {
        if(vd_remove_set.find(*n_it) == vd_remove_set.end())
        {
            //log_node_info(*n_it, true, false, false);
            vd_t vd_new = local_condense(*n_it, vd_remove_set);
            if (vd_new != NULL)
            {
                w_list.push_back(vd_new);
                //log_node_info(vd_new, true, false, false);
            }
            else
            {
                //shc_log_info(shc_logname, "no condense \n");
            }
        }
    }
    num_condensed_rm += vd_remove_set.size();


    //shc_log_info(shc_logname, "remove nodes associated with local condense\n");
    for(std::set<vd_t>::iterator it=vd_remove_set.begin();
                            it != vd_remove_set.end(); it++)
    {
        //shc_log_info(shc_logname, "in %d, out %d\n", boost::in_degree(*it, graph),
        assert(boost::in_degree(*it, graph)==0 &&
                boost::out_degree(*it, graph)==0);
        boost::remove_vertex(*it,graph);
        assert(1==exist_vd.erase(*it));
    }
}

bool Sequence_graph_handler::
is_read_bridge_node(read_num_t info_read_id, read_length_t info_start, vd_t vd)
{
    Read_acc acc = coll_read_list.get_read(info_read_id);
    if (info_start==0)
        return false;
    if (acc.len <= info_start + graph[vd].seq_len())
    {
        return false;
    }
    bool flag = memcmp(&graph[vd].seq.at(0),
            acc.read_ptr+info_start, graph[vd].seq_len())==0;
    return flag;
}

vd_t Sequence_graph_handler::
create_left_node(ed_t ed)
{
    //shc_log_info(shc_logname, "Start create_left_node\n");
    vd_t vd_source = boost::source(ed, graph);
    vd_t vd_target = boost::target(ed, graph);
    bundled_node_p & source_node = graph[vd_source];
    bundled_node_p & curr_node = graph[vd_target];

    vd_t u = boost::add_vertex(graph);
    exist_vd.insert(u);
    if(vd_source != vd_target) // not self node
    {
        aer_t new_aer = boost::add_edge(vd_source, u, graph);
        graph[new_aer.first] = bundled_edge_p(graph[ed].weight + 1, graph[ed].count);
    }

    //shc_log_info(shc_logname, "source %s\n", source_node.seq.c_str());
    //shc_log_info(shc_logname, "target %s\n", curr_node.seq.c_str());

    int left_index = source_node.seq.size() - 1 - graph[ed].weight;
    //shc_log_info(shc_logname, "left index %d\n", left_index);

    double count = (source_node.count + curr_node.count)/2;
    std::string seq = source_node.seq.at(left_index) + curr_node.seq;
    //shc_log_info(shc_logname, "seq %s\n", seq.c_str());
    graph[u] = bundled_node_p(seq, glob_node_id++);

    //debug
    //if(vd_source == vd_target)
    //{
    //    shc_log_info(shc_logname, "SELF LOOP\n");
    //    shc_log_info(shc_logname, "%d %s\n", left_index, seq.c_str());
    //    shc_log_info(shc_logname, "source %s\n", source_node.seq.c_str());
    //    shc_log_info(shc_logname, "target %s\n", curr_node.seq.c_str());
    //}
    //std::vector<Bdg_read_info>::iterator it;
    //for(it=curr_node.reads_info.begin(); it!=curr_node.reads_info.end(); ++it)
    //{
    //    if(is_read_bridge_node(it->read_id, it->start-1, u))
    //    {
    //        node_link_read(u, it->read_id, it->start-1);
    //    }
    //}
    //shc_log_info(shc_logname, "Finish create_left_node\n");
    return u;
}

vd_t Sequence_graph_handler::
create_right_node(ed_t ed)
{
    vd_t vd_source = boost::source(ed, graph);
    vd_t vd_target = boost::target(ed, graph);
    bundled_node_p & curr_node = graph[vd_source];
    bundled_node_p & target_node = graph[vd_target];

    vd_t w = boost::add_vertex(graph);
    exist_vd.insert(w);
    if(vd_source != vd_target) //not self node
    {
        aer_t new_aer = boost::add_edge(w, boost::target(ed, graph), graph);
        graph[new_aer.first] = bundled_edge_p(graph[ed].weight + 1, graph[ed].count);
    }

    std::string seq = curr_node.seq + target_node.seq.at(graph[ed].weight);
    graph[w] = bundled_node_p(seq, glob_node_id++);

    //std::vector<Bdg_read_info>::iterator it;
    //for(it=curr_node.reads_info.begin(); it!=curr_node.reads_info.end(); ++it)
    //{
    //    if(is_read_bridge_node(it->read_id, it->start+1, w))
    //    {
    //        node_link_read(w, it->read_id, it->start+1);
    //    }
    //}
    return w;
}

int Sequence_graph_handler::find_known_path(int max_hop)
{
    start_timer(&timer);
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start find known path\n");
#endif

    typedef tsl::hopscotch_map<uint64_t, std::vector<Node_index>,
                   hash_u64, equ64> Kmer_ni_Map;
    typedef tsl::hopscotch_map<uint64_t, std::vector<Node_index>,
                    hash_u64, equ64>::iterator Kmer_ni_Map_iterator;
    Kmer_ni_Map kmer_ni_map;

    std::string kmer_base;
    kmer_base.resize(kmer_length);

    // construct kmer to node-index multimap, it is needed
    // since we perform multibridging, there can be multiple vertex that
    // that contain the start of a read if it is not a bridging read
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Building kmer->vd start_index multimap\n");
#endif
    uint64_t byte;
    vip_t vip = boost::vertices(graph);
    for(vi_t vi=vip.first; vi!=vip.second; ++vi)
    {
        bundled_node_p & curr_node = graph[*vi];
        std::string & node_seq = curr_node.seq;
        for(int i=0; i<curr_node.seq_len()-kmer_length+1; ++i)
        {
            encode_kmer(&node_seq.at(i), &byte, kmer_length);
            Kmer_ni_Map_iterator it = kmer_ni_map.find(byte);
            if(it != kmer_ni_map.end())
            {
                (it.value()).emplace_back(*vi, i);
            }
            else
            {
                std::vector<Node_index> vec_seq(1, Node_index(*vi, i));
                vec_seq.reserve(KNOWN_PATH_VEC_INIT_SIZE);
                kmer_ni_map.insert(std::make_pair(byte, vec_seq));
            }
        }
    }

    int num_known_path = 0;
    int num_known_edge = 0;
    int num_checked_read = 0;
    int num_searched_reads = 0;
    uint64_t start_byte;
    uint64_t end_byte;
    //iterate through all read to get all known path
    //shc_log_info(shc_logname, "Iterating reads to find path\n"); << std::endl;
    //std::ofstream writer(setting.local_files.output_path+"/all_reads.log");

    /*
    coll_read_list.clear();

    std::string read_path_single_prefix =
        "/data1/bowen/Shannon_D_seq/output_SE_WW_R/back/back_all_reads/components_reads/comp";
    size_t file_size =
        get_total_num_reads(read_path_single_prefix);
    //std::cout << "comp " << comp_i << "single file size " << file_size << std::endl;
    coll_read_list.s_reads.setup(file_size);
    Kmer_Node_map kmer_node_map;
    int i=0;
    std::string file_path = read_path_single_prefix + std::to_string(curr_comp) + "_" + std::to_string(i);
    while(exist_path(file_path))
    {
        load_all_single_read(file_path, kmer_node_map);
        i++;
        file_path = read_path_single_prefix + "_" + std::to_string(i);
    }
    coll_read_list.declare_read_finish();
    */

    //std::cout << "num reads " << coll_read_list.get_num_reads() << std::endl;
    for(read_num_t i=0; i<coll_read_list.get_num_reads(); i++)
    {
        //shc_log_info(shc_logname, "start new read\n");
        num_checked_read++;
        Read_acc acc = coll_read_list.get_read(i);
        //std::string a_read(acc.read_ptr, acc.len);
        //writer <<a_read << std::endl;
        encode_kmer(acc.read_ptr, &start_byte, kmer_length);
        encode_kmer(acc.read_ptr+acc.len-kmer_length, &end_byte, kmer_length);

#ifdef LOG_SEQ_GRAPH
        char start_back[33];
        char end_back[33];
        start_back[kmer_length] = '\0';
        end_back[kmer_length] = '\0';
        decode_kmer( start_back, &start_byte, kmer_length);
        decode_kmer( end_back, &end_byte, kmer_length);
        std::string read_seq(acc.read_ptr, acc.len);
        shc_log_info(shc_logname, "Read  %s\n", read_seq.c_str());
        shc_log_info(shc_logname, "start %s\n", start_back);
        shc_log_info(shc_logname, "end   %s\n", end_back);
#endif

        if(kmer_ni_map.find(end_byte) == kmer_ni_map.end())
            continue;
        Kmer_ni_Map_iterator start_it = kmer_ni_map.find(start_byte);
        if(start_it == kmer_ni_map.end())
            continue;


        int curr_hop = max_hop;

        std::vector<Node_index> & nodes = start_it.value();

        //shc_log_info(shc_logname, "search matched nodes\n");
        num_searched_reads++;
        for(int j=0; j<nodes.size(); j++)
        {
            vd_t vd = nodes[j].vd;
            read_length_t val_node_index = nodes[j].start;
            bundled_node_p & start_node = graph[vd];
            char * val_node_ptr = &start_node.seq.at(val_node_index);
            read_length_t val_node_len = start_node.seq_len() - val_node_index;

            // if the first node matches the start of the read
            if(is_partial_match(acc.read_ptr, acc.len,
                                val_node_ptr, val_node_len))
            {
                std::vector<vd_t> path;
                if(search_seq(acc.read_ptr,
                              acc.len,
                              vd,
                              val_node_index, curr_hop, path))
                {
                    // there is a full cover path
                    //shc_log_info(shc_logname, "find matched path\n");
                    //read_known_path_map[i] = path;

                    for(uint64_t i=0; i<path.size()-1; i++)
                    {
                        aer_t aer = boost::edge(path[i], path[i+1], graph);
                        assert(aer.second);
                        std::map<ed_t, uint64_t>::iterator e_it =
                                                  known_edges.find(aer.first);
                        if(e_it == known_edges.end())
                        {
                            //shc_log_info(shc_logname, "read count %d\n", acc.read_count);
                            known_edges.insert(std::make_pair(aer.first, acc.read_count));
                            num_known_edge++ ;
                        }
                        else
                        {
                            e_it->second += acc.read_count;
                        }
                    }

                    if(path.size() > 2)
                    {
                        known_path_set.insert(path);
                    }
                }
            }
        }
    }
    //std::cout << "num_checked_read " << num_checked_read << std::endl;
    //std::cout << "num_searched_reads " << num_searched_reads << std::endl;

    num_known_path = known_path_set.size();
    //std::cout << "find " << num_known_path << " known paths\n";
    //std::cout << "find " << num_known_edge << " known edges\n";
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Finish find known path, Finding %d full matching path\n", num_known_path);
#endif
    //std::cout << "finding " << num_known_path << " known path finish" << std::endl;
#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "finding " << num_known_path << " known path finish, ";
    stop_timer(&timer);
#endif

    return num_known_path;
    //bebug
    /*
    for(std::set<std::vector<vd_t> >::iterator it=known_path_set.begin();
                                it!=known_path_set.end(); it++)
    {
        const std::vector<vd_t> & a_path = *it;
        shc_log_info(shc_logname, "new path size %d\n", a_path.size());
        for(std::vector<vd_t>::const_iterator p_it=a_path.begin();
                            p_it!=a_path.end(); p_it++)
        {
            shc_log_info(shc_logname, "%s\n", graph[*p_it].seq.c_str());
        }
    }
    */
}

void Sequence_graph_handler::log_classify_edge_types()
{
    std::ofstream writer(setting.local_files.output_path + "/prune_edge");

    int case_1 = 0;
    int case_2 = 0;
    int case_3 = 0;
    int total_0_edge = 0;
    int total_edge = 0;
    int total_case_1 = 0;
    int total_case_2 = 0;
    int total_case_3 = 0;
    eip_t eip = boost::edges(graph);
    for(ei_t ei=eip.first; ei!=eip.second; ei++)
    {
        vd_t vd_1 = boost::source(*ei, graph);
        vd_t vd_2 = boost::target(*ei, graph);
        int in_degree_v1 = boost::in_degree(vd_1, graph);
        int out_degree_v2 = boost::out_degree(vd_2, graph);
        int out_degree_v1 = boost::out_degree(vd_1, graph);
        int in_degree_v2 = boost::in_degree(vd_2, graph);
        if (graph[*ei].count == 0)
        {
            total_0_edge++;
            if(in_degree_v1==0 || out_degree_v2==0)
                case_1 ++;
            else if (out_degree_v1==1 || in_degree_v2 == 1)
                case_2 ++;
            else if (out_degree_v1>=2 and in_degree_v2 >=2)
                case_3 ++;
            else
            {
                std::cout << "other cases "
                          <<  " in_degree_v1 " << in_degree_v1
                          <<  " out_degree_v1 " << out_degree_v1
                          <<  " in_degree_v2 " << in_degree_v2
                          <<  " out_degree_v2 " << out_degree_v2 << std::endl;
                exit(0);
            }
        }

        total_edge++;
        if(in_degree_v1==0 || out_degree_v2==0)
            total_case_1 ++;
        else if (out_degree_v1==1 || in_degree_v2 == 1)
            total_case_2 ++;
        else if (out_degree_v1>=2 and in_degree_v2 >=2)
            total_case_3 ++;
        else
        {
            std::cout << "other cases "
                      <<  " in_degree_v1 " << in_degree_v1
                      <<  " out_degree_v1 " << out_degree_v1
                      <<  " in_degree_v2 " << in_degree_v2
                      <<  " out_degree_v2 " << out_degree_v2 << std::endl;
            exit(0);
        }
    }
    writer << "total_0_edge\t" << total_0_edge << "\t" << total_edge << std::endl;
    writer << "case_1\t" << case_1 << "\t" << total_case_1 << std::endl;
    writer << "case_2\t" << case_2 << "\t" << total_case_2 << std::endl;
    writer << "case_3\t" << case_3 << "\t" << total_case_3 << std::endl;
}

void Sequence_graph_handler::log_classify_node_types(std::ofstream & writer, int id)
{
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        std::string & node_seq = graph[vd].seq;
        int out_degree = boost::out_degree(vd, graph);
        int in_degree = boost::in_degree(vd, graph);

        writer << id << "\t" << node_seq << "\t" << in_degree << "\t"
               << out_degree;

        in_eip_t in_eip = boost::in_edges(vd,graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
           vd_t vd_source = boost::source(*in_ei, graph);
           writer << "\t" << graph[vd_source].seq << "\t" << graph[*in_ei].weight
                  << "\t" << graph[*in_ei].count;
        }

        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            writer << "\t" << graph[vd_target].seq << "\t" << graph[*out_ei].weight
                   << "\t" << graph[*out_ei].count;
        }
        writer << std::endl;
    }
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
//#ifdef LOG_SEQ_GRAPH
//        if(val_node_len!=0)
//            path_seq.insert(path_seq.end(), curr_node.seq.begin()+val_start, curr_node.seq.end());
//#endif
        if(curr_seq_len <= val_node_len)
        {
//#ifdef LOG_SEQ_GRAPH
//            shc_log_info(shc_logname, "success %u, %u\n", curr_seq_len, val_node_len);
//            path_seq.resize(seq_full_len);
//            shc_log_info(shc_logname, "         %s\n", path_seq.c_str());
//#endif
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
    return memcmp(ptr1, ptr2, std::min(len1, len2))==0; //is thread sage by POSIX
}

void Sequence_graph_handler::break_self_loops()
{
    start_timer(&timer);
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "start break self loops\n");
#endif

    std::vector<out_ei_t> self_loop_edge;

    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            if(vd_target == vd)  // a self loop
            {
                self_loop_edge.push_back(out_ei);
                //bundled_node_p & node = graph[vd];
                //bundled_edge_p & edge_p = graph[*out_ei];
                //node.seq = node.seq + node.seq.substr(edge_p.weight);
                if(boost::out_degree(vd, graph)==1)
                {
                    shc_log_info(shc_logname, "have a isolated self loops\n");
                }
                //break; //since for each node, it only contains one self node
            }
        }
    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "have %d self loops\n", self_loop_edge.size());
#endif
    for(std::vector<out_ei_t>::iterator l_it=self_loop_edge.begin();
                    l_it!=self_loop_edge.end(); l_it++)
    {
        boost::remove_edge(*l_it, graph);
    }

#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "break all self loop finish, ";
    stop_timer(&timer);
#endif
}

void Sequence_graph_handler::py_break_cycles()
{
    std::set<vd_t> no_cycle;

    std::deque<vd_t> cycle;
    while(py_find_cycle(no_cycle, cycle))
    {
        simple_break_cycle(cycle);
        cycle.clear();
    }
    condense_graph();
}

// return true when detect a cycle
bool Sequence_graph_handler::py_find_cycle(std::set<vd_t> & no_cycle,
            std::deque<vd_t> & cycle)
{
    shc_log_info(shc_logname, "start find cycle\n");
    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        std::set<vd_t>::iterator nc_it = no_cycle.find(vd);
        if(nc_it == no_cycle.end())
        {
            std::vector<vd_t> traversed;
            py_reachable_cycle(vd, no_cycle, traversed, cycle);
            if (cycle.size() > 0)
                return true;
        }
    }
    return false;
}

void Sequence_graph_handler::py_reachable_cycle( vd_t vd,
            std::set<vd_t> & no_cycle, std::vector<vd_t> & traversed, std::deque<vd_t> & cycle)
{
    traversed.push_back(vd);
    shc_log_info(shc_logname, "curr vd %u\n", graph[vd].node_id);
    for(int i=0; i<traversed.size(); i++)
    {
        shc_log_info(shc_logname, "traversed %d %u\n",i, graph[traversed[i]].node_id);
    }
    out_eip_t out_eip = boost::out_edges(vd,graph);
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        std::vector<vd_t>::iterator it =
                std::find(traversed.begin(), traversed.end(), vd_target);
        if(it != traversed.end()) // contain cycle
        {
            cycle.clear();
            for(std::vector<vd_t>::iterator tr_it=it; tr_it!=traversed.end(); tr_it++)
            {
                cycle.push_back(*tr_it);
            }
            cycle.push_back(vd_target);
            shc_log_info(shc_logname, "detect a cycle\n");
            for(int i=0; i<cycle.size(); i++)
            {
                shc_log_info(shc_logname, "cycle %d %u\n",i, graph[cycle[i]].node_id);
            }
            return;
        }
        if (no_cycle.find(vd) != no_cycle.end())
            continue;
        py_reachable_cycle( vd_target, no_cycle, traversed, cycle);
        if(cycle.size() > 0)
            return;
    }
    no_cycle.insert(vd);
    return;
}

void Sequence_graph_handler::break_all_cycles()
{
    start_timer(&timer);
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start break all cycle\n");
#endif
    std::set<vd_t> acyclic_node_set;
    std::deque<vd_t> cycle_path;
    int i=0;
    int j=0;
    while (find_cycle(acyclic_node_set, cycle_path))
    {
        //shc_log_info(shc_logname, "Start break cycle with size %d\n", cycle_path.size());

        simple_break_cycle(cycle_path);
        //if(!break_cycle(cycle_path))
        //{
        //    simple_break_cycle(cycle_path);
        //}
        //shc_log_info(shc_logname, "start break path\n");
        //
        cycle_path.clear();
        j++;
        //shc_log_info(shc_logname, "after clear path\n");
        //std::cout << "finish break " << (++i) << " cycles" << std::endl;
    }

    //std::cout << "Finish break " << j << " cycle\n";
    //log_graph_to_file(graph_writer , 200);
    condense_graph();
    //condense_graph();
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Finish break %d cycle\n", j);
#endif
#ifdef SHOW_SEQ_GRAPH_PROGRESS
    std::cout << "break all cycles finish, ";
    stop_timer(&timer);
#endif

}

/**
 * Delete first node
 * @param cycle_path
 */
void Sequence_graph_handler::simple_break_cycle(std::deque<vd_t> & cycle_path)
{
    vd_t vd_remove = cycle_path[1]; //according to python implementation, take 1
    //std::cout << graph[vd_remove].seq << std::endl;
    //shc_log_info(shc_logname, "rm %u\n", graph[vd_remove].node_id);
    //boost::remove_edge(cycle_path[0], cycle_path[1], graph);
    boost::clear_vertex(vd_remove, graph);
    assert(1==exist_vd.erase(vd_remove));
    boost::remove_vertex(vd_remove, graph);
}

bool Sequence_graph_handler::
find_cycle(std::set<vd_t> & acyclic_node_set, std::deque<vd_t> & cycle_path)
{
    //shc_log_info(shc_logname, "start find_cycle\n");
    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        std::set<vd_t>::iterator an_it = acyclic_node_set.find(vd);
        if(an_it == acyclic_node_set.end()) //not in the set
        {
            if(is_node_inside_cycle(vd, acyclic_node_set, cycle_path))
            {
                //shc_log_info(shc_logname, "finish find_cycle\n");
                for(int i=0; i<cycle_path.size(); i++)
                    cycle_writer << graph[cycle_path[i]].seq << "\t";
                cycle_writer << std::endl;
                return true;
            }
        }
    }
    //shc_log_info(shc_logname, "finish find_cycle\n");
    return false;
}

std::deque<vd_t>::iterator Sequence_graph_handler::
thread_safe_find(std::deque<vd_t> & cycle_path, vd_t vd_target)
{
    for(std::deque<vd_t>::iterator it= cycle_path.begin();
                                   it!=cycle_path.end(); it++)
    {
        if(*it==vd_target)
            return it;
    }
    return cycle_path.end();
}

/**
 * Search cycle in depth first fashion
 * return true if cycle is detected, false otherwise
 */
bool Sequence_graph_handler::
is_node_inside_cycle(vd_t vd, std::set<vd_t> & acyclic_node_set,
                                            std::deque<vd_t> & cycle_path)
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
        cycle_index_it = thread_safe_find(cycle_path, vd_target);
        //cycle_index_it = std::find(cycle_path.begin(), cycle_path.end(), vd_target);
        //if find same vertex in the middle
        if(cycle_index_it  != cycle_path.end())
        {
            int mid_index = std::distance(cycle_path.begin(), cycle_index_it);
            //shc_log_info(shc_logname, "%d find in middle %d, node id %u\n",
            //                curr_comp, mid_index, graph[*cycle_index_it].node_id);

            cycle_path.push_back(vd_target); //  path now contain cycle
            //shc_log_info(shc_logname, "cycle node %s\n", graph[vd_target].seq.c_str());
            if(cycle_index_it != cycle_path.begin() )
            {
                /*
                shc_log_info(shc_logname, "%d cycle path ", curr_comp);
                for(std::deque<vd_t>::iterator cp_it= cycle_path.begin();
                                               cp_it!=cycle_path.end(); cp_it++)
                {
                    vd_t vd = *cp_it;
                    assert(exist_vd.find(vd)!=exist_vd.end());
                    bundled_node_p & curr_node = graph[vd];
                    info_log_info(shc_logname, "%u ",
                            curr_node.node_id);
                }
                */
                //shc_log_info(shc_logname, "\n");
                for(int i=0; i<mid_index; i++)
                {
                    //vd_t erase_vd = cycle_path.front();
                    //shc_log_info(shc_logname, "before erase %d, node id %u\n",
                    //        i, graph[erase_vd].node_id);
                    cycle_path.pop_front();
                    //shc_log_info(shc_logname, "after erase one\n");
                }

                //shc_log_info(shc_logname, "after erase all\n");
            }
            else
            {
                //shc_log_warning("cycle path contain two identical nodes \n");
            }
            // now the cycle_path variable contains just the cycle

            //assert((cycle_path.back()) == (cycle_path.front()));
            return true;
        }
        if(acyclic_node_set.find(vd_target) != acyclic_node_set.end())
            continue;

        if(is_node_inside_cycle(vd_target, acyclic_node_set, cycle_path))
            return true;
    }
    acyclic_node_set.insert(vd);
    cycle_path.pop_back();
    //cycle_path.clear();
    //shc_log_info(shc_logname, "QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n");
    return false;
}

/**
 * Used for checking final answer if the entire graph is acyclic
 * @return false when not acyclic
 */
bool Sequence_graph_handler::acyclic_check()
{
    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        S_info & dfs_info = graph[vd].s_info;
        dfs_info.color = WHITE;
        dfs_info.parent = static_cast<void *>(NULL);
    }

    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        S_info & dfs_info = graph[vd].s_info;
        if(dfs_info.color == WHITE)
        {
            if(!DFS_visit(vd))
            {
#ifdef LOG_SEQ_GRAPH
                shc_log_info(shc_logname, "cycle node vd %d\n", graph[vd].node_id);
#endif
                return false;
            }
        }
    }
    return true;
}

// return false when it is not acyclic
// implement from <<introduction to algorithm>>
bool Sequence_graph_handler::DFS_visit(vd_t vd)
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
            v_dfs_info.parent = (void *)vd;
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

// return true if cycle is properly handled
bool Sequence_graph_handler::break_cycle(std::deque<vd_t> cycle_path)
{
    // find anything to break
    assert(cycle_path.size()>=2);
    /*
    std::cout << "cycle_path size " << cycle_path.size()  << std::endl;
    for(std::deque<vd_t>::iterator c_it=cycle_path.begin();
                                    c_it!=cycle_path.end(); ++c_it)
    {
        std::cout <<  graph[*c_it].node_id << " " ;
    }
    std::cout <<  std::endl;
    */
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
                prev_node.push_back(cycle_path.back());
            else
                prev_node.push_back(*(c_it-1));
            if(c_it==cycle_path.end()-1)
                next_node.push_back(cycle_path.front());
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
                    prev_node.push_back(cycle_path.back());
                else
                    prev_node.push_back(*(c_it-1));
                if(c_it==cycle_path.end()-1)
                    next_node.push_back(cycle_path.front());
                else
                    next_node.push_back(*(c_it+1));
            }
        }
    }
    if(xnodes.size()==0)
    {
        shc_log_info(shc_logname, "no enough xnodes in cycle, unable to break\n");
        /*
        std::cout << "cycle_path size " << cycle_path.size()  << std::endl;
        for(std::deque<vd_t>::iterator c_it=cycle_path.begin();
                                        c_it!=cycle_path.end(); ++c_it)
        {
            std::cout <<  graph[*c_it].node_id << " " ;
        }
        std::cout << std::endl;
        */
        return false;
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
                o_in_copy_count_sum += graph[*in_ei].count;
            else
                c_in_copy_count = graph[*in_ei].count;
        }

        int  o_out_copy_count_sum=0, c_out_copy_count=0;
        out_eip_t out_eip = boost::out_edges(vd_curr, graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            if(vd_target != vd_next)
                o_out_copy_count_sum += graph[*out_ei].count;
            else
                c_out_copy_count = graph[*out_ei].count;
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
    return true;
}

void Sequence_graph_handler::
break_cycle_splitting(vd_t vd_split, vd_t vd_c_in, vd_t vd_c_out,
                                    int sum_copy_count, int sum_c_in_and_o_out)
{
    // actual breaking
    //shc_log_info(shc_logname, "breaking cycle node %s\n", graph[vd_split].seq.c_str());
    vd_t vd_clone = boost::add_vertex(graph);
    exist_vd.insert(vd_clone);
    double proportion = 0.5;
    if ( sum_copy_count != 0 )
    {
        proportion = static_cast<double>(sum_c_in_and_o_out) /
                     static_cast<double>(sum_copy_count);
    }
    double clone_count = (proportion*graph[vd_split].prevalence);

    graph[vd_clone] = bundled_node_p(graph[vd_split].seq, glob_node_id++);
    graph[vd_clone].reads_info = graph[vd_split].reads_info;
    graph[vd_split].prevalence = ((1-proportion)*graph[vd_split].prevalence);

    aer_t aer_cycle_in = boost::edge(vd_c_in, vd_split, graph);
    assert(aer_cycle_in.second);
    aer_t aer_cycle_out = boost::edge(vd_split, vd_c_out, graph);
    assert(aer_cycle_out.second);

    aer_t aer_new = boost::add_edge(vd_c_in, vd_clone, graph);
    graph[aer_new.first] = bundled_edge_p(graph[aer_cycle_in.first].weight,
                                          graph[aer_cycle_in.first].count);
    boost::remove_edge(aer_cycle_in.first, graph);

    out_eip_t out_eip = boost::out_edges(vd_split, graph);
    std::vector<ed_t> edge_to_remove;
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        if(vd_target != vd_c_out)
        {
            aer_new = boost::add_edge(vd_clone, vd_target, graph);
            graph[aer_new.first] = bundled_edge_p(graph[aer_cycle_out.first].weight,
                                        graph[aer_cycle_out.first].count);
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

// check all in edges from vd1, and link them to vd
void Sequence_graph_handler::link_in_edges(vd_t vd1, vd_t vd)
{
    aer_t aer;
    in_eip_t vd1_in_eip = boost::in_edges(vd1, graph);
    for(in_ei_t it = vd1_in_eip.first; it!=vd1_in_eip.second; it++)
    {
        aer = boost::add_edge(boost::source(*it, graph), vd, graph);
        graph[aer.first] = bundled_edge_p(graph[*it].weight);
    }
    //boost::clear_in_edges(vd1,graph);
}

void Sequence_graph_handler::link_out_edges(vd_t vd2, vd_t vd)
{
    aer_t aer;
    out_eip_t vd2_out_eip = boost::out_edges(vd2, graph);
    for(out_ei_t it = vd2_out_eip.first; it!=vd2_out_eip.second; it++)
    {
        aer = boost::add_edge(vd, boost::target(*it, graph), graph);
        graph[aer.first] = bundled_edge_p(graph[*it].weight);
    }
    //boost::clear_out_edges(vd2,graph);
}

void Sequence_graph_handler::seq_graph_output_dir_setup_helper(std::string & dir)
{
    boost::filesystem::path dir_boost_path(dir);
    if(boost::filesystem::exists(dir))
        empty_directory(dir_boost_path, setting.local_files.output_path);
    else
        add_directory(dir_boost_path);
}

void Sequence_graph_handler::output_components(std::string & node_dir,
                                std::string & edge_dir, std::string & path_dir)
{
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "Start output_components\n");
#endif
    //create dir

    seq_graph_output_dir_setup_helper(node_dir);
    seq_graph_output_dir_setup_helper(edge_dir);
    seq_graph_output_dir_setup_helper(path_dir);
    std::string single_node_path(node_dir + "/single_node_" + std::to_string(curr_comp));
    //std::cout << "write to " << single_node_path << std::endl;
    std::ofstream sinlge_node_writer(single_node_path);

    typedef std::multimap<vd_t, std::vector<vd_t> > Paths_by_start_MMap;
    typedef std::multimap<vd_t, std::vector<vd_t> >::iterator
                                             Paths_by_start_MMap_iterator;
    typedef std::pair<Paths_by_start_MMap_iterator, Paths_by_start_MMap_iterator>
                                        Paths_by_start_MMap_iterator_pair;
    Paths_by_start_MMap paths_by_start_mmap;
    for(std::set<std::vector<vd_t> >::iterator it= known_path_set.begin();
                                               it!=known_path_set.end(); ++it)
    {
        paths_by_start_mmap.insert(std::make_pair(it->at(0),*it));
    }

    std::set<vd_t> black_nodes;
    uint64_t i =0;
    uint64_t single_seq_i = 0;

    uint64_t num_edge_written = 0;
    uint64_t num_edge_skipped = 0;
    //for all nodes
    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd_curr = *it;
        if(black_nodes.find(vd_curr) == black_nodes.end())
        {
            std::set<vd_t> nodes;
            std::vector<vd_t> sorted_nodes;
            std::set<ed_t> edges;

            add_component(vd_curr, nodes, edges);

#ifdef LOG_SEQ_GRAPH
            shc_log_info(shc_logname, "nodes has size %d\n", nodes.size());
#endif
            bool is_DAC = variant_topological_sort(nodes, sorted_nodes);
            assert(is_DAC);
            assert(sorted_nodes.size()==nodes.size());

#ifdef LOG_SEQ_GRAPH
            shc_log_info(shc_logname, "after sort, nodes has size %d\n", sorted_nodes.size());
#endif
            //shc_log_info(shc_logname, "after sort, nodes has size %d\n", sorted_nodes.size());
            if(sorted_nodes.size()==1 || edges.size()==0)
            {
                //shc_log_info(shc_logname, "ID: %d %s\n", graph[sorted_nodes.at(0)].node_id,
                //                        graph[sorted_nodes.at(0)].seq.c_str());

                black_nodes.insert(vd_curr);
                //std::cout << "graph[vd_curr].seq_len() " << graph[vd_curr].seq_len() << std::endl;
                if(graph[vd_curr].seq_len() >
                                setting.output_seq_min_len)
                {
                    sinlge_node_writer << ">s_" << curr_comp << "_"
                                 << (single_seq_i++) << std::endl;
                    sinlge_node_writer << graph[vd_curr].seq << std::endl;
                }
                continue;
            }
            std::ofstream node_file(node_dir + "/node" + std::to_string(i));
            std::ofstream edge_file(edge_dir + "/edge" + std::to_string(i));
            std::ofstream path_file(path_dir + "/path" + std::to_string(i));

            // output nodes
            node_file << "ID\tBases\tCount" << std::endl;
            for(std::vector<vd_t>::iterator it= sorted_nodes.begin();
                                            it!=sorted_nodes.end(); ++it)
            {
                node_file << graph[*it].node_id << "\t" + graph[*it].to_string() << std::endl;
                black_nodes.insert(*it);
            }

            path_file << "VD1\tVD2 ..." << std::endl;
            // output path

            for(std::vector<vd_t>::iterator it= sorted_nodes.begin();
                                            it!=sorted_nodes.end(); ++it)
            {
                Paths_by_start_MMap_iterator_pair ps_it_pair =
                                        paths_by_start_mmap.equal_range(*it);

                for(Paths_by_start_MMap_iterator ps_it =ps_it_pair.first;
                                           ps_it!=ps_it_pair.second; ps_it++)
                {
                    std::string path_string;
                    path_to_string(ps_it->second, path_string);
                    path_file << path_string << std::endl;
                }
            }

            // output edge
            edge_file << "VD\tVD\tweight\tcount" << std::endl;
            for(std::set<ed_t>::iterator it=edges.begin(); it!=edges.end(); ++it)
            {
                if(graph[*it].count > 0)
                {
                    std::string temp = edge_to_string(*it);
                    //shc_log_info(shc_logname, "test %s\n", temp.c_str());
                    edge_file << edge_to_string(*it) << std::endl;
                    num_edge_written ++;
                }
                else
                {
                    num_edge_skipped++;
                }
            }
            i++;

            node_file.close();
            edge_file.close();
            path_file.close();
        }
    }
    //std::cout << "num_edge_skipped " << num_edge_skipped << std::endl;
    //std::cout << "num_edge_written " << num_edge_written << std::endl;
    num_single_seq = single_seq_i;
    sinlge_node_writer.close();
    if(!is_run_single_component)
    {
        std::string cmd("cp " + single_node_path + " " +
                    setting.local_files.single_node_dir + "/single_node_" + std::to_string(curr_comp));
        //std::cout << cmd << std::endl;
        run_command(cmd, false);
    }
}

std::string Sequence_graph_handler::
edge_to_string(ed_t ed)
{
    vd_t vd_source = boost::source(ed, graph);
    vd_t vd_target = boost::target(ed, graph);
    return std::to_string(graph[vd_source].node_id) + "\t" +
           std::to_string(graph[vd_target].node_id) + "\t" +
           std::to_string(graph[ed].weight) + "\t" +
           std::to_string(graph[ed].count);

}

void Sequence_graph_handler::
path_to_string(std::vector<vd_t> & path, std::string & out_string)
{
    out_string.clear();
    for(std::vector<vd_t>::iterator it=path.begin(); it!=path.end(); ++it)
    {
        out_string += (std::to_string(graph[*it].node_id) + '\t');
    }
}

// return false when the graph is not a DAC
bool Sequence_graph_handler::
variant_topological_sort(std::set<vd_t> & nodes, std::vector<vd_t> & sorted_nodes)
{
    std::set<vd_t> added;
    std::stack<vd_t> fringe;

    for(std::set<vd_t>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        vd_t vd_curr = *it;
        if(boost::in_degree(vd_curr, graph) == 0 )
            fringe.push(vd_curr);
    }

    while (fringe.size()>0)
    {
        vd_t vd = fringe.top();
        fringe.pop();
        if (added.find(vd) != added.end())
            continue;
        added.insert(vd);
        sorted_nodes.push_back(vd);
        out_eip_t out_eip = boost::out_edges(vd, graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            if(is_all_predecessors_included(vd_target, added))
                fringe.push(vd_target);
        }
    }
    if (added.size() != nodes.size())
        return false;
    else
        return true;
}

bool Sequence_graph_handler::
is_all_predecessors_included(vd_t vd_root, std::set<vd_t> & added_nodes)
{
    in_eip_t in_eip = boost::in_edges(vd_root,graph);
    for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
    {
        vd_t vd_source = boost::source(*in_ei, graph);
        if(added_nodes.find(vd_source) == added_nodes.end())
            return false;
    }
    return true;
}


void Sequence_graph_handler::
add_component(vd_t vd_root, std::set<vd_t> & nodes, std::set<ed_t> & edges)
{
    std::stack<vd_t> grey_nodes;
    grey_nodes.push(vd_root);
    while (grey_nodes.size()>0)
    {
        vd_t vd = grey_nodes.top();
        grey_nodes.pop();
        if (nodes.find(vd) != nodes.end())
            continue;

        nodes.insert(vd);

        out_eip_t out_eip = boost::out_edges(vd, graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            edges.insert(*out_ei);
            grey_nodes.push(boost::target(*out_ei, graph));
        }

        in_eip_t in_eip = boost::in_edges(vd, graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
            //edges.insert(*in_ei);  so that there is no cycle
            grey_nodes.push(boost::source(*in_ei, graph));
        }
    }
}

// return true if it is suspicious
bool Sequence_graph_handler::is_suspicious(vd_t vd)
{
    bundled_node_p & curr_node = graph[vd];
    int inD = boost::in_degree(vd, graph);
    int outD = boost::out_degree(vd, graph);
    //std::cout << "size_threshold " << size_threshold << std::endl;
    //std::cout << "prevalence_threshold " << prevalence_threshold << std::endl;
    if((inD==0 || outD==0) && curr_node.seq_len()<=size_threshold)
        return true;

    if(average_prevalence(vd) >= prevalence_threshold)
        return false;

    if (inD==0 || outD==0 )
    {
        return true;
    }

    //std::cout << "reach is suspicuous last section" << std::endl;
    return true;
}

// return false for no suspicious
bool Sequence_graph_handler::remove_suspicious_nodes(uint64_t & node_removed)
{
    //shc_log_info(shc_logname, "start remove_suspicious_nodes\n");
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    std::vector<Node_Count_Pair> suspicous_nodes;
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        if(is_suspicious(*it))
        {
            //shc_log_info(shc_logname, "Suspicious\n");
            suspicous_nodes.emplace_back(*it, average_prevalence(*it));
        }
    }

    if(suspicous_nodes.size()==0)
        return false;

    Node_count_sorter sorter;
    std::sort(suspicous_nodes.begin(), suspicous_nodes.end(), sorter);


    std::set<vd_t> vd_remove_set;

    for(std::vector<Node_Count_Pair>::iterator it=suspicous_nodes.begin();
                                               it!=suspicous_nodes.end(); it++ )
    {
        std::set<vd_t> adjacent_vd;
        vd_t vd = it->first;
        if(boost::in_degree(vd,graph)==0 && boost::out_degree(vd,graph)==0)
        {
            vd_remove_set.insert(vd);
            continue;
        }

        in_eip_t in_eip = boost::in_edges(vd,graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
            vd_t vd_source = boost::source(*in_ei, graph);
            adjacent_vd.insert(vd_source);
        }

        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            adjacent_vd.insert(vd_target);
        }
        boost::clear_vertex(vd, graph);
        vd_remove_set.insert(vd);
        //boost::remove_vertex(vd, graph);
        for(std::set<vd_t>::iterator adj_it= adjacent_vd.begin();
                                        adj_it!=adjacent_vd.end(); ++adj_it)
        {
            local_condense(*adj_it, vd_remove_set);
        }
    }
    node_removed = vd_remove_set.size();
    remove_set_nodes(vd_remove_set);
    return true;
}

void Sequence_graph_handler::remove_all_suspicious_nodes()
{
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "remove all suspicous\n");
#endif
    uint64_t node_removed = 0;
    uint64_t total_node_removed = 0;
    while(remove_suspicious_nodes(node_removed))
    {
        total_node_removed += node_removed;
    }
    //std::cout << total_node_removed << " suspicous node removed\n" << std::endl;

#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "%u suspicous node removed\n", total_node_removed);
    shc_log_info(shc_logname, "finish remove all suspicous\n");
#endif
}

void Sequence_graph_handler::collapse_all()
{
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "collapse all nodes\n");
#endif
    uint64_t node_removed = 0;
    uint64_t total_node_removed = 0;
    while(collapse_nodes(node_removed))
    {
        total_node_removed += node_removed;
    }
#ifdef LOG_SEQ_GRAPH
    shc_log_info(shc_logname, "%u node collapsed\n", total_node_removed);
    shc_log_info(shc_logname, "finish all nodes\n");
#endif
}

bool Sequence_graph_handler::collapse_nodes(uint64_t & node_removed)
{
    bool is_collapse = false;
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    std::set<vd_t> vd_remove;
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        if(collapse_successor(*it, vd_remove))
        {
            is_collapse = true;
        }
    }

    node_removed = vd_remove.size();
    // remove nodes
    if (is_collapse)
        remove_set_nodes(vd_remove);
    return is_collapse;
}

void Sequence_graph_handler::remove_set_nodes(std::set<vd_t> & vd_remove)
{
    for(std::set<vd_t>::iterator it= vd_remove.begin();
                                 it!=vd_remove.end(); ++it)
    {
        std::vector<Terminal_node_info> & term_node_info = graph[*it].term_nodes_info;
        for(std::vector<Terminal_node_info>::iterator t_it=term_node_info.begin();
                    t_it!=term_node_info.end(); t_it++ )
        {
            p1_end[t_it->align_read_id] = NULL;
            p2_front[t_it->align_read_id] = NULL;
        }
        boost::clear_vertex(*it, graph);
        boost::remove_vertex(*it, graph);
        assert(1==exist_vd.erase(*it));
    }
}

bool Sequence_graph_handler::collapse_successor(vd_t vd, std::set<vd_t> & vd_remove)
{
    out_eip_t out_eip = boost::out_edges(vd, graph);
    std::vector<vd_t> successor;
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        successor.push_back(vd_target);
    }

    for(int i=0; i<successor.size(); i++)
    {
        for(int j=i+1; j<successor.size(); j++)
        {
            if(similar(successor[i], successor[j]))
            {
                vd_remove.insert(collapse(successor[i], successor[j]));
                return true;
            }
        }
    }
    return false;
}

bool Sequence_graph_handler::similar(vd_t vd1, vd_t vd2)
{
    bundled_node_p & node_1 = graph[vd1];
    bundled_node_p & node_2 = graph[vd2];
    if(node_1.seq_len() != node_2.seq_len())
        return false;

    double num_diff = 0;
    int base_len = node_1.seq_len();
    for(int i=0; i<base_len; i++)
    {
        if(node_1.seq[i] != node_2.seq[i])
            num_diff += 1.0f;
    }

    if(num_diff/static_cast<double>(base_len) >= hamming_frac)
        return false;

    std::set<vd_t> successor_vd1;
    std::set<vd_t> successor_vd2;
    get_successors(vd1, successor_vd1);
    get_successors(vd2, successor_vd2);
    if(successor_vd1 != successor_vd2)
        return false;
    std::set<vd_t> predecessor_vd1;
    std::set<vd_t> predecessor_vd2;
    get_predecessors(vd1, predecessor_vd1);
    get_predecessors(vd2, predecessor_vd2);
    if(predecessor_vd1 != predecessor_vd2)
        return false;

    return true;
}

void Sequence_graph_handler::erase_terminal_node_info(vd_t vd)
{
    std::vector<Terminal_node_info> & term_nodes_info = graph[vd].term_nodes_info;
    for(std::vector<Terminal_node_info>::iterator it=term_nodes_info.begin();
                                        it!=term_nodes_info.end(); it++)
    {
        p1_end[it->align_read_id] = NULL;
        p2_front[it->align_read_id] = NULL;
    }
}


// return nodes to be removed
vd_t Sequence_graph_handler::collapse(vd_t vd_1, vd_t vd_2)
{
    bundled_node_p & node_1 = graph[vd_1];
    bundled_node_p & node_2 = graph[vd_2];
    //std::cout << "prevalence 1 " << node_1.prevalence << " "<< node_1.seq <<std::endl;
    //std::cout << "prevalence 2 " << node_2.prevalence << " "<< node_2.seq <<std::endl;
    if(node_1.prevalence > node_2.prevalence)
    {
        node_1.prevalence +=  node_2.prevalence;
        erase_terminal_node_info(vd_2);
        boost::clear_vertex(vd_2, graph);

        //std::cout << "collapse 2 " << node_2.prevalence << std::endl;
        return vd_2;
    }
    else if(node_1.prevalence < node_2.prevalence)
    {
        node_2.prevalence += node_1.prevalence;
        erase_terminal_node_info(vd_1);
        boost::clear_vertex(vd_1, graph);
        //std::cout << "collapse 1 " << node_1.prevalence << std::endl;
        return vd_1;
    }
    else
    {
        if(node_1.seq > node_2.seq)
        {
            node_2.prevalence += node_1.prevalence;
            erase_terminal_node_info(vd_1);
            boost::clear_vertex(vd_1, graph);
            //std::cout << "collapse " << node_1.seq << std::endl;
            return vd_1;
        }
        else
        {
            node_1.prevalence +=  node_2.prevalence;
            erase_terminal_node_info(vd_2);
            boost::clear_vertex(vd_2, graph);
            //std::cout << "collapse " << node_2.seq << std::endl;
            return vd_2;
        }


    }
}

void Sequence_graph_handler::get_successors(vd_t vd, std::set<vd_t> & successors)
{
    out_eip_t out_eip = boost::out_edges(vd,graph);
    for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
    {
        vd_t vd_target = boost::target(*out_ei, graph);
        successors.insert(vd_target);
    }
}

void Sequence_graph_handler::
get_predecessors(vd_t vd, std::set<vd_t> & predecessors)
{
    in_eip_t in_eip = boost::in_edges(vd, graph);
    for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
    {
        vd_t vd_source = boost::source(*in_ei, graph);
        predecessors.insert(vd_source);
    }
}

void Sequence_graph_handler::find_copy_count()
{
    int check_num_known_edges = 0;
    vip_t vip = boost::vertices(graph);
    std::set<ed_t> remove_edges;
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        //uint64_t total_node = 0;
        vd_t vd_source = *it;
        bundled_node_p & node_source = graph[vd_source];
        out_eip_t out_eip = boost::out_edges(vd_source, graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            double ec = 0;
            bundled_edge_p & edge_p = graph[*out_ei];
            std::map<ed_t, uint64_t>::iterator e_it = known_edges.find(*out_ei);

            vd_t vd_1 = boost::source(*out_ei, graph);
            vd_t vd_2 = boost::target(*out_ei, graph);
            int in_degree_v1 = boost::in_degree(vd_1, graph);
            int out_degree_v2 = boost::out_degree(vd_2, graph);
            int out_degree_v1 = boost::out_degree(vd_1, graph);
            int in_degree_v2 = boost::in_degree(vd_2, graph);

            //if(edge_p.count<=1  &&
            //  (in_degree_v1 <=2 || out_degree_v2<=2 || out_degree_v1<=2 || in_degree_v2<=2))
            //{
            //    ec = 0;
            //    remove_edges.insert(*out_ei);
            //}

            if(e_it == known_edges.end() ) // no such edge
            {
                /*


                if(in_degree_v1 != 0 && out_degree_v2 != 0)
                {
                    if (out_degree_v1==1 || in_degree_v2 == 1)
                    {
                        ec = 0;
                        remove_edges.insert(*out_ei);
                    }
                }
                continue; // no update edge count
                */
                ec = 0;
                remove_edges.insert(*out_ei);

            }
            else
            {
                ec = e_it->second;
                assert(e_it->second != 0);
                check_num_known_edges++;
            }

            //total_node += ec
            edge_p.count = ec/std::max(
                    static_cast<double>(coll_read_list.L- edge_p.weight-1),
                    static_cast<double>(1));

        }
    }
    for(std::set<ed_t>::iterator it=remove_edges.begin();
        it != remove_edges.end(); it++)
    {
        boost::remove_edge(*it, graph);
    }

    //std::cout << "check num known edges " << check_num_known_edges << std::endl;
}

void Sequence_graph_handler::find_approximate_copy_count()
{
    //shc_log_info(shc_logname, "find_approximate_copy_count\n");
    vip_t vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        bundled_node_p & node = graph[vd];
        node.norm = node.seq_len() - kmer_length + 1;
        node.copy_count = double(node.prevalence)/node.norm;
    }

    vip = boost::vertices(graph);
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd_source = *it;
        bundled_node_p & node_source = graph[vd_source];
        out_eip_t out_eip = boost::out_edges(vd_source, graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            bundled_node_p & node_target = graph[vd_target];
            bundled_edge_p edge_p = graph[*out_ei];
            kmer_count_node_t norm = std::max(
                    static_cast<kmer_count_node_t>(coll_read_list.L-edge_p.weight-1),
                    static_cast<kmer_count_node_t>(0));

            kmer_count_node_t node_counts =
                                    node_source.copy_count * node_source.norm +
                                    node_target.copy_count * node_target.norm;
            if (norm==0)
                edge_p.count = 0;
            else
                edge_p.count = 0.5 *static_cast<double>(node_counts) /
                                    static_cast<double>(norm);
        }
    }
}

// return true if all nodes in this read exists in the filtered kmer dict
vd_t Sequence_graph_handler::
update_edge_count_with_read(std::string & base, Kmer_Node_map & kmer_node_map)
{
    Kmer_Node_map_iterator it_end = kmer_node_map.end();
    Kmer_Node_map_iterator it_prev, it_next;
    int i=0;
    bool find_next = false;
    Kmer_Node_map_iterator it_curr = kmer_node_map.find(base.substr(i, kmer_length));
    //shc_log_info(shc_logname, "%d:%s %s\n", i, base.substr(i, kmer_length).c_str(),((it_curr!=it_end)?("Yes"):("No")));
    while(i<base.size()-kmer_length)
    {
        ++i;
        //shc_log_info(shc_logname, "current i: %d\n", i);
        if(it_curr!=it_end)
        {
            it_next = kmer_node_map.find(base.substr(i, kmer_length));
            //shc_log_info(shc_logname, "%d: %s %s\n", i, base.substr(i, kmer_length).c_str(),((it_curr!=it_end)?("Yes"):("No")));
            if(it_next != it_end)
            {
                aer_t aer = boost::edge(it_curr->second, it_next->second, graph);
                graph[aer.first].count += 1;
            }
            it_curr = it_next;
        }
        else
        {
            it_curr = kmer_node_map.find(base.substr(i, kmer_length));
            //shc_log_info(shc_logname, "%d:%s %s\n", i, base.substr(i, kmer_length).c_str(),((it_curr!=it_end)?("Yes"):("No")));
        }
    }
}


void Sequence_graph_handler::remove_nodes_edges_if_not_cover_by_reads()
{
    /*
    if(setting.filter_single)
    {
        shc_log_info(shc_logname, "total num nodes %d\n", get_num_nodes());
        vip_t vip = boost::vertices(graph);
        //walk through the graph list, and update if
        std::set<vd_t> remove_vds;
        for(vi_t it=vip.first; it!=vip.second; it++)
        {
            vd_t & vd = *it;
            bundled_node_p & node = graph[vd];
            if( node.reads_info.size() == 0 &&
                node.seq_len() < min_read_len-2)
            {
                shc_log_info(shc_logname, "min len %u\n", min_read_len);
                boost::clear_vertex(vd, graph);
                remove_vds.insert(vd);
            }
        }
        for(std::set<vd_t>::iterator it=remove_vds.begin(); it!=remove_vds.end(); it++)
        {
            boost::remove_vertex(*it, graph);
        }
    }
    */
}



void Sequence_graph_handler::log_all_node(bool is_log_node, bool is_log_read)
{
    std::ofstream debug_file(setting.local_files.output_path+"/all_node.log");
    shc_log_info(shc_logname, "log all nodes, number of vertices is %d, number of edges is %d\n",
                boost::num_vertices(graph), boost::num_edges(graph));
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        log_node_info(*it, is_log_node, is_log_read, false);
        //if(graph[*it].seq_len() > 24)
        //debug_file << (graph[*it].seq) << " " << boost::in_degree(*it, graph) << " "  <<boost::out_degree(*it, graph)<< std::endl;
            //debug_file << (graph[*it].seq) << " " << graph[*it].reads_info.size() << std::endl;
        /*
        in_eip_t in_eip = boost::in_edges(*it,graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
            vd_t vd_source = boost::source(*in_ei, graph);
            debug_file << (graph[*it].seq) << " " << graph[vd_source].seq << std::endl;
        }
        out_eip_t out_eip = boost::out_edges(*it,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            debug_file << (graph[*it].seq) << " " << graph[vd_target].seq<< std::endl;
        }
        */
    }
    debug_file.close();
}

void Sequence_graph_handler::log_xnode_set(bool is_log_node, bool is_log_read)
{
    std::ofstream debug_file(setting.local_files.output_path+"/xnode.log");
    shc_log_info(shc_logname, "log all xnodes\n");
    for(std::set<vd_t>::iterator it=xnode_set.begin(); it!=xnode_set.end(); it++)
    {
        log_node_info(*it, is_log_node, is_log_read, false);
        debug_file << (graph[*it].seq) << std::endl;
    }
    debug_file.close();
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

void Sequence_graph_handler::log_node_info(vd_t vd, bool is_log_nodes ,
                                        bool is_log_read_info, bool is_log_bfs)
{
    bundled_node_p & curr_node = graph[vd];
    assert(!curr_node.seq.empty());
    assert(boost::in_degree(vd,graph)<=4 && boost::out_degree(vd,graph)<=4);
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
            shc_log_info(shc_logname, "%d %d %d %d %s\n", boost::in_degree(vd_source, graph),
                    boost::out_degree(vd_source, graph), graph[vd_source].node_id,
                    graph[*in_ei].weight, graph[vd_source].seq.c_str());
        }
        shc_log_info(shc_logname, "Out nodes:\n");
        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            shc_log_info(shc_logname, "%d %d %d %d %s\n", boost::in_degree(vd_target, graph),
                    boost::out_degree(vd_target, graph) , graph[vd_target].node_id,
                    graph[*out_ei].weight ,graph[vd_target].seq.c_str());
        }
    }

    shc_log_info(shc_logname, "%d Read cover it\n", curr_node.reads_info.size());
    if(is_log_read_info)
    {
        std::vector<Bdg_read_info> & reads_info = curr_node.reads_info;
        for(std::vector<Bdg_read_info>::iterator it=reads_info.begin();
                it!=reads_info.end(); it++)
        {
            //shc_log_info(shc_logname, "total_num read %u\n", coll_read_list.get_num_read());

            shc_log_info(shc_logname, "read id %u, start %u\n", it->read_id, it->start);
            log_read_seq(coll_read_list.get_read(it->read_id));
        }
    }
    assert(boost::in_degree(vd,graph)<=4 && boost::out_degree(vd,graph)<=4);

    if(is_log_bfs)
    {
        log_node_bfs_info(vd, true, true);
    }
    shc_log_info(shc_logname,"\n");
}

void Sequence_graph_handler::log_all_edge(bool is_log_weight, bool is_log_count)
{
    shc_log_info(shc_logname, "log all edges, number of edges is %d\n",
                                                       boost::num_edges(graph));
    eip_t eip = boost::edges(graph);
    for(ei_t ei=eip.first; ei!=eip.second; ei++)
    {
        if(is_log_weight)
            shc_log_info(shc_logname, "weight %d\n", graph[*ei].weight);
        if(is_log_count)
            shc_log_info(shc_logname, "count 4.2%\n", graph[*ei].count);

    }
}

void Sequence_graph_handler::
log_node_bfs_info(vd_t vd, bool is_show_color, bool is_shown_depth)
{
    S_info & bfs_info = graph[vd].s_info;
    if(is_show_color)
    {
        if(bfs_info.color == WHITE)
            shc_log_info(shc_logname, "color WHITE\n");
        else if(bfs_info.color == GRAY)
            shc_log_info(shc_logname, "color GRAY\n");
        else if(bfs_info.color == BLACK)
            shc_log_info(shc_logname, "color BLACK\n");
        else
        {
            shc_log_info(shc_logname, "color Unknown\n");
            exit(1);
        }
    }
    if(is_shown_depth)
    {
        shc_log_info(shc_logname, "depth %d\n", bfs_info.d);
    }
}

void Sequence_graph_handler::log_edge_count()
{
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        shc_log_info(shc_logname, "******************************************\n");
        shc_log_info(shc_logname, "VD: %d, Seq: %s\n", graph[vd].node_id, graph[vd].seq.c_str());
        in_eip_t in_eip = boost::in_edges(vd,graph);
        for(in_ei_t in_ei=in_eip.first; in_ei!=in_eip.second; in_ei++)
        {
            vd_t vd_source = boost::source(*in_ei, graph);
            shc_log_info(shc_logname, "IN  VD: %6d, weight %6d, count %4.2f\n",
                graph[vd_source].node_id, graph[*in_ei].weight, graph[*in_ei].count);
        }
        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            shc_log_info(shc_logname, "OUT VD: %6d, weight %6d, count %4.2f\n",
                graph[vd_target].node_id, graph[*out_ei].weight, graph[*out_ei].count);
        }
    }
}

void Sequence_graph_handler::log_term_array(bool is_show_seq)
{
    shc_log_info(shc_logname, "\n");
    shc_log_info(shc_logname, "\t\t LOG TERM NODE\n");
    assert(p1_end.size()== p2_front.size());
    for(uint64_t i=0; i<p1_end.size(); i++)
    {
        vd_t vd1 = p1_end[i];
        vd_t vd2 = p2_front[i];
        shc_log_info(shc_logname, "\n");
        if(vd1 != NULL && vd2!=NULL)
        {
            std::set<vd_t>::iterator vd1_it = exist_vd.find(vd1);
            std::set<vd_t>::iterator vd2_it = exist_vd.find(vd2);

            if(vd1_it==exist_vd.end())
            {
                shc_log_info(shc_logname, "align index %u: p1 term node do not exist\n", i);
                continue;
            }
            if(vd2_it==exist_vd.end())
            {
                shc_log_info(shc_logname, "align index %u: p2 term node do not exist\n", i);
                continue;
            }

            shc_log_info(shc_logname, "align index %u: node id: %u %u\n", i,
                                       graph[vd1].node_id, graph[vd2].node_id);

            if(is_show_seq)
            {
                shc_log_info(shc_logname, "p1 seq: %s\n", graph[vd1].seq.c_str());
                shc_log_info(shc_logname, "p2 seq: %s\n", graph[vd2].seq.c_str());
            }
        }
        else if(vd1 == NULL && vd2==NULL)
        {
            shc_log_info(shc_logname, "%u: No terminal node exists\n", i);
        }
        else
        {
            shc_log_error("only one node exists\n");
            exit(0);
        }

    }
    shc_log_info(shc_logname, "\t\t END LOG TERM NODE\n\n");
}

int Sequence_graph_handler::get_num_isolated_nodes()
{
    vip_t vip = boost::vertices(graph);
    int num_isolated_nodes = 0;
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        if(boost::in_degree(vd, graph) == 0 &&
           boost::out_degree(vd, graph) == 0)
        {
            num_isolated_nodes++;
        }
    }
    return num_isolated_nodes;
}

void Sequence_graph_handler::log_graph_to_file(std::ofstream & writer, int id)
{
    vip_t vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        std::string & node_seq = graph[vd].seq;
        int out_degree = boost::out_degree(vd, graph);
        int in_degree = boost::in_degree(vd, graph);
        out_eip_t out_eip = boost::out_edges(vd,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            writer << id << "\t" << node_seq << "\t" << graph[vd_target].seq << "\t"
                   << in_degree << "\t" << out_degree << "\t"
                   << graph[*out_ei].weight << std::endl;
        }

        if( out_degree==0 && boost::in_degree(vd, graph)==0)
            writer << id << "\t" << node_seq << std::endl;

    }
}

void Sequence_graph_handler::log_all_nodes_reads_to_file(std::ofstream & writer)
{
    vip_t vip = boost::vertices(graph);
    writer << "algo iter 0" << std::endl;
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd = *it;
        bundled_node_p & curr_node = graph[vd];
        writer << "seq " << curr_node.seq.c_str()  << std::endl;
        writer << "in  degree " << boost::in_degree(vd, graph) << std::endl;
        writer << "out degree " << boost::out_degree(vd, graph) << std::endl;
        writer << "num read " << curr_node.reads_info.size() << std::endl;
        for(std::vector<Bdg_read_info>::iterator ri_it=curr_node.reads_info.begin();
                                          ri_it!=curr_node.reads_info.end(); ++ri_it)
        {
            // if the same read is processed, ignore the rest
            Read_acc acc = coll_read_list.get_read(ri_it->read_id);
            std::string read(acc.read_ptr, acc.len);
            writer << "read " <<   read << std::endl;
        }
    }
}
