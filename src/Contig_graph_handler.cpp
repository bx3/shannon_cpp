#include "Contig_graph_handler.h"

uint8_t get_num_bit[ 256 ] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

size_t contig_curr_proc = 0;
size_t contig_prev_proc = 0;

uint64_t total_length = 0;
uint64_t total_query=0 ;
uint64_t total_update=0 ;
uint64_t conn_size = 0;
uint64_t total_bit_set_query = 0;

struct Block_timer cgh_timer;
struct Block_timer local_timer;
struct Block_timer metis_timer;

struct Block_timer get_set_timer;


void Contig_graph_handler::run_contig_graph_handler()
{
    // partition contig
    struct Block_timer cgh_main_timer;
    start_timer(&cgh_main_timer);

    group_components();

    //shc_log_info(shc_logname, "log complex comp num edges\n");
    //for(int i=0; i<comps_num_edges.size();i++ )
    //{
    //    shc_log_info(shc_logname, "comp %d has %d edges\n", i, comps_num_edges[i]);
    //}

    shc_log_info(shc_logname, "finish group contig\n");
    log_stop_timer(&cgh_main_timer);

    //std::cout << "finish group contig, ";
    //stop_timer(&cgh_main_timer);

    dump_component_array(setting.local_files.output_comp_path);

    add_or_overwrite_directory(lf.output_components_read_dir, lf.output_path);
    add_or_overwrite_directory(lf.output_components_read_features_dir, lf.output_path);



    dump_filtered_contig();

    fasta_dumper.setup_dump_files(component_size);

    //std::ofstream ana_file_writer(lf.output_path + "/sampling.log");
    comp_total_kmer_count.assign(curr_component_num, std::vector<uint64_t>());
    comp_num_reads.assign(curr_component_num, 0);
    comp_num_kmers.assign(curr_component_num, 0);
    read_sampler.setup(curr_component_num, setting.random_seed);


    // assigning reads
    if(setting.has_single)
    {
        if(read_sampler.is_store_all_reads_and_features)
        {
            single_read_features_dumper.setup_dump_files(component_size, 0);
            single_read_features_dumper.setup_disk_dumper(lf.output_components_read_features_dir,
                                                          lf.comp_read_features_prefix);
        }

        start_timer(&cgh_main_timer);
        //assign_reads_to_components(setting.local_files.input_read_path);
        assign_reads_to_components_mmap(setting.local_files.input_read_path);
        if(read_sampler.is_store_all_reads_and_features)
        {
            single_read_features_dumper.write_mem_to_disk();
            single_read_features_dumper.reset();
        }
        shc_log_info(shc_logname, "finish assigning single read\n");
        log_stop_timer(&cgh_main_timer);
#ifdef SHOW_PROGRESS
#endif
    }
    if(setting.has_pair)
    {
        start_timer(&cgh_timer);
        if(read_sampler.is_store_all_reads_and_features)
        {
            paired_read_features_dumper.setup_dump_files(component_size, 0);
            paired_read_features_dumper.setup_disk_dumper(lf.output_components_read_features_dir,
                                                          lf.comp_pair_read_features_prefix);
        }
        assign_paired_read_to_components_mmap(
                            setting.local_files.input_read_path_1,
                            setting.local_files.input_read_path_2);
        if(read_sampler.is_store_all_reads_and_features)
        {
            paired_read_features_dumper.write_mem_to_disk();
            paired_read_features_dumper.reset();
        }
        //assign_paired_read_to_components(
        //                    setting.local_files.input_read_path_1,
        //                    setting.local_files.input_read_path_2);
        shc_log_info(shc_logname, "finish assigning paired read\n");
        log_stop_timer(&cgh_timer);
#ifdef SHOW_PROGRESS
#endif
    }
    fasta_dumper.finalize_dump_files();


    // dumping kmers
    start_timer(&cgh_timer);
    kmer_dumper.setup_dump_files(component_size, NUM_OPEN_FILE);
    assign_kmer_to_components();
    shc_log_info(shc_logname, "finish assigning kmer\n");
    log_stop_timer(&cgh_timer);
#ifdef SHOW_PROGRESS
#endif

    //for(uint64_t i=0; i<curr_component_num ; i++)
    //{
    //    for(uint64_t j=0; j<comp_total_kmer_count[i].size(); j++ )
    //    {
    //        ana_file_writer << comp_total_kmer_count[i][j] << "\t";
    //    }
    //    ana_file_writer << comp_num_reads[i]
    //                    << "\t" << comp_num_kmers[i] << std::endl;
    //}
    //ana_file_writer.close();
}

void Contig_graph_handler::dump_filtered_contig()
{
    std::ofstream writer(setting.local_files.filtered_contig);
    uint64_t write_index = 0;
    for(uint64_t i=0; i<component_array.size(); i++)
    {
        comp_num_t comp_j = component_array[i];
        if(complex_comp_indicator[comp_j])
        {
            writer << ">Contig_" << write_index << std::endl;
            char * base_start;
            uint64_t len;
            ch->get_contig(i, base_start, len);
            std::string contig(base_start, len);
            writer << contig << std::endl;
            write_index ++;
        }
    }
}

//return true when graph is needed, false if it is ok
bool Contig_graph_handler::
assign_comp_without_graph()
{
    while(!contig_stack.empty())
    {
        contig_num_t i = contig_stack.back();
        conn_contig_set.set_conn_contig(i);
        contig_stack.pop_back();
        char * base_start;
        uint64_t len;
        ch->get_contig(i, base_start, len);
        //shc_log_info(shc_logname, "new contig\n");
        total_length += len;
        //vd_t curr_vd;
        //assert(contig_vertex_map.get_contig_vd(i, curr_vd));

        for(Contig_handler::size_type j=0; j<len-kmer_length+1;   j++)
        {
            //if(j%setting.contig_graph_setup.down_sample == 0)
            //{
            uint64_t byte;
            encode_kmer(base_start+j, &byte, kmer_length);
            Kmer_counter_map_iterator it = kh->kmer_counter.find(byte);
            //if(it == kh->kmer_counter.end())
            //{
            //    char temp[33];
            //    memcpy(temp, base_start+j, kmer_length);
            //    temp[kmer_length] = '\0';
            //    std::cout << "cannot find kmer " << temp << std::endl;
            //    exit(0);
            //}

            if(j==0 || j==len-kmer_length)
            {
                if(get_num_bit[it->second.info] == 1)
                    continue;
            }
            else
            {
                if(get_num_bit[it->second.info] == 2)
                    continue;
            }

            char next_letter, before_letter;
            if(j==0)
            {
                before_letter = ANY_LETTER;
            }
            else
            {
                before_letter = *(base_start+j -1);
            }

            if(j!=len-kmer_length)
            {
                next_letter = *(base_start+j+kmer_length);
            }
            else
            {
                next_letter = ANY_LETTER;
            }
            get_connected_contig_index_without_graph(it, contig_stack,
                            before_letter, next_letter);
        }
        if(conn_contig_set.num_contig > metis_setup.partition_size)
        {
            for(contig_num_t i=0; i<conn_contig_set.num_contig; i++ )
            {

                //shc_log_warning("partition_size %u\n",(metis_setup.partition_size));
                //shc_log_warning("num contig %u\n",(conn_contig_set.num_contig));
                //shc_log_warning("%u\n",(conn_contig_set.connected_contig)[i]);
                explorable_contig_set.erase_explored(
                            (conn_contig_set.connected_contig)[i]);
            }
            conn_contig_set.reset();

            return true; // graph is needed
        }
    }
    conn_size += conn_contig_set.num_contig;
    //std::cout << "conn_contig_set size " << conn_contig_set.size() << std::endl;
#ifdef LOG_CONTIG_GRAPH
    shc_log_info(shc_logname, "no metis\n");
#endif
    if(!is_set_collect_comp_num)
    {
        collect_comp_num = curr_component_num;
        std::string simple_comp(SIMPLE_COMPONENT);
        component_type.push_back(simple_comp);
        complex_comp_indicator.push_back(false);
        accum_collect_contig_num = 1;
        curr_component_num++;
        is_set_collect_comp_num = true;
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "start comp %d\n", curr_component_num);
#endif
    }

    //update contig to component info
    for(contig_num_t i=0; i<conn_contig_set.num_contig; i++ )
    {
        //shc_log_warning("%u\n",(conn_contig_set.connected_contig)[i]);
        component_array[(conn_contig_set.connected_contig)[i]] = collect_comp_num;
        accum_collect_contig_num++;
    }

    // start with a new component number
    if (accum_collect_contig_num >= non_partition_size)
    {
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "component %u has %d no-metis contig\n",
                        collect_comp_num, accum_collect_contig_num);
#endif
        is_set_collect_comp_num = false;
#ifdef LOG_CONTIG_GRAPH
    shc_log_info(shc_logname, "start a new collect component\n");
#endif
    }
    conn_contig_set.reset();
    return false;
}


void Contig_graph_handler::group_components()

{
    shc_log_info(shc_logname, "Start constructing contig graph\n");
    char kmer_array[33];
    kmer_array[kmer_length] = '\0';
    uint64_t byte = 0;
    Block_timer progress_timer;
    start_timer(&progress_timer);


    Block_timer test_timer;
    uint64_t num_graph = 0;
    uint64_t num_non_graph = 0;
    uint64_t num_conn = 0;

    start_timer(&cgh_timer);

    {
        std::string message = "Building contig graph using " + std::to_string(ch->num_contig)
                            + " contigs\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t progress_step = ch->num_contig/100;

        while(!explorable_contig_set.is_explore_all())
        {
#ifdef LOG_CONTIG_GRAPH
            shc_log_info(shc_logname, "\t\t\t\t\t#processing component %u\n",
                                                                curr_component_num);
#endif
#ifdef SHOW_PROGRESS
            if((contig_curr_proc-contig_prev_proc) > progress_step)
            {
                int percentage = (100 * contig_curr_proc)/(ch->num_contig);
                progress.write(static_cast<double>(contig_curr_proc) /
                                    static_cast<double>(ch->num_contig));
                //std::cout << "[" << percentage << "%] " <<  contig_curr_proc
                     //<< " contig out of " <<   ch->num_contig
                     //<< " is processed, contig len processed " << total_length
                     //<< ", total_query " << total_query
                     //<< ", num graph " << num_graph
                     //<< ", num non graph " << num_non_graph
                     //<< ", conn_size "  << conn_size
                     //<< ", total_bit_set_query " << total_bit_set_query
                     //<< ", total_update " << total_update
                     //<< std::endl;

                shc_log_info(shc_logname, "[ %d %], %d  contig out of %d "
                            "is processed, contig len processed %d, total_query %d, "
                            "num graph %u, num non graph %u, num conn_size %u, "
                            "total_bit_set_query %u, total_update %u \n",
                            percentage, contig_curr_proc, ch->num_contig, total_length,
                            total_query, num_graph, num_non_graph, conn_size,
                            total_bit_set_query, total_update);
                conn_size = 0;
                num_graph = 0;
                num_non_graph = 0;
                total_bit_set_query = 0;
                contig_prev_proc = contig_curr_proc;
                total_update = 0;

                log_stop_timer(&progress_timer);
                start_timer(&progress_timer);
                total_query=0;
                total_length = 0;
            }
#endif
            //accurate_start_timer(&get_set_timer);
            contig_num_t root_contig = explorable_contig_set.get_set_next_explorable_contig();
            //accurate_accumulate_timer(&get_set_timer);
            contig_stack.clear();
            contig_stack.push_back(root_contig);

            if(assign_comp_without_graph()) //
            {
                num_graph++;
                contig_stack.clear();
                contig_stack.push_back(root_contig);
                assign_comp_with_graph();
            }
            else
            {
                num_non_graph++;
            }
            num_conn ++;
            contig_curr_proc = explorable_contig_set.num_explored_contig;
        }
    }

    count_component_size();

    //std::cout << "Finish finding " << curr_component_num << " component, ";
    //stop_timer(&cgh_timer);
    shc_log_info(shc_logname, "Finish constructing contig graph\n");
}

void Contig_graph_handler::assign_comp_with_graph()
{
    // add that node, so each graph has at least one node
    Contig_vertex_map contig_vertex_map(ch->num_contig);
    graph_t graph;
    vd_t vd = boost::add_vertex(graph);
    contig_vertex_map.set_contig_vd(contig_stack.back(), vd);
    graph[vd] = bundled_contig_index(contig_stack.back());
    //accurate_start_timer(&test_timer);
    //explore that contig and find new contig
    while(!contig_stack.empty())
    {
        contig_num_t i = contig_stack.back();
        contig_stack.pop_back();
        //curr_contig = i;
        char * base_start;
        uint64_t len;
        ch->get_contig(i, base_start, len);

        //shc_log_info(shc_logname, "new contig\n");
        total_length += len;

        vd_t curr_vd;
        assert(contig_vertex_map.get_contig_vd(i, curr_vd));

        for(Contig_handler::size_type j=0;
            j<len-kmer_length+1;   j++)
        {
            //if(j%setting.contig_graph_setup.down_sample == 0)
            //{
            //conn_point = j;
            uint64_t byte;
            encode_kmer(base_start+j, &byte, kmer_length);
            Kmer_counter_map_iterator it = kh->kmer_counter.find(byte);

            if(j==0 || j==len-kmer_length)
            {
                if(get_num_bit[it->second.info] == 1)
                    continue;
            }
            else
            {
                if(get_num_bit[it->second.info] == 2)
                    continue;
            }

            char next_letter, before_letter;
            if(j==0)
                before_letter = ANY_LETTER;
            else
                before_letter = *(base_start+j -1);
            if(j!=ch->contig_len_list.at(i)-kmer_length)
                next_letter = *(base_start+j+kmer_length);
            else
                next_letter = ANY_LETTER;
            get_connected_contig_index(it, graph, contig_vertex_map, contig_stack,
                            curr_vd, before_letter, next_letter);
        }
    }
    break_and_keep_component(graph);
    graph.clear();
}

void Contig_graph_handler::count_component_size()
{
    component_size.resize(curr_component_num, 0);

    for(int i=0; i<component_array.size(); i++)
    {
        Contig_handler::size_type contig_len = ch->contig_len_list.at(i);
        comp_num_t comp_i = component_array[i];
        comp_num_t comp_i_aux = component_array_aux[i];
        if(comp_i < curr_component_num) //is valid components
        {
            component_size[comp_i] += contig_len;
        }
        if(comp_i_aux < curr_component_num) //is valid components
        {
            component_size[comp_i_aux] += contig_len;
        }
    }
    //for(std::vector<uint64_t>::iterator it=component_size.begin();
    //                it!=component_size.end(); it++)
    //{
    //    shc_log_info(shc_logname, "%u\n", *it);
    //}

    if(setting.has_single)
    {
        num_read_files += curr_component_num;
        num_kmer_files += curr_component_num;
    }
    if(setting.has_pair)
    {
        num_read_files += curr_component_num*2;
        num_kmer_files += curr_component_num;
    }
}

void Contig_graph_handler::break_and_keep_component(graph_t & graph)
{
    typedef std::pair<vi_t, vi_t> vip_t;
    idx_t num_vertices = boost::num_vertices(graph);
    //shc_log_info(shc_logname, "num vertices %u, partition_size %u\n",num_vertices , metis_setup.partition_size);
    if( num_vertices <= metis_setup.partition_size)
    {
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "no metis\n");
#endif
        shc_log_error("should not enter here\n");
        exit(1);

        if(!is_set_collect_comp_num)
        {
            collect_comp_num = curr_component_num;
            std::string simple_comp(SIMPLE_COMPONENT);
            component_type.push_back(simple_comp);
            complex_comp_indicator.push_back(false);
            accum_collect_contig_num = 1;
            curr_component_num++;
            is_set_collect_comp_num = true;
#ifdef LOG_CONTIG_GRAPH
            shc_log_info(shc_logname, "start comp %d\n", curr_component_num);
#endif
        }

        std::pair<vi_t, vi_t> vip = boost::vertices(graph);
        //update contig to component info
        for(vi_t it=vip.first; it!=vip.second; it++)
        {
            component_array[graph[*it].contig_index] = collect_comp_num;
            accum_collect_contig_num++;
        }
        // start with a new component number
        if (accum_collect_contig_num >= non_partition_size)
        {
            shc_log_info(shc_logname, "component %u has %d no-metis contig\n",
                            collect_comp_num, accum_collect_contig_num);

            shc_log_info(shc_logname, "num_edges %d, num nodes %d\n",
                     boost::num_edges(graph), boost::num_vertices(graph));
            is_set_collect_comp_num = false;
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "start a new collect component\n");
#endif
        }
    }
    else        //call metis
    {
//#ifdef LOG_METIS
        shc_log_info(shc_logname, "Starting metis \n");
//#endif
        shc_log_info(shc_logname, "num_edges %d, num nodes %d\n",
                     boost::num_edges(graph), boost::num_vertices(graph));

        create_metis_array_input(graph);
        metis_setup.num_vertices = num_vertices;

        idx_t partition = static_cast<idx_t>(
                std::ceil(static_cast<float>(num_vertices)/
                          static_cast<float>(metis_setup.partition_size)));    //note
        //std::cout << "num vertices " << num_vertices << std::endl;
        //std::cout << "partition size " << metis_setup.partition_size << std::endl;

        shc_log_info(shc_logname, "partition to %d\n", partition);

        metis_setup.num_partition = std::min(static_cast<idx_t>(100), partition);
        //shc_log_info(shc_logname, "metis to asked to break %d partition\n", metis_setup.num_partition);
        idx_t ufactor = static_cast<idx_t>(METIS_IMBALANCE_PRECISION
                                                   * (metis_setup.overload-1));
        METIS_SetDefaultOptions(metis_setup.options);
        metis_setup.options[METIS_OPTION_UFACTOR] = ufactor;
        metis_setup.options[METIS_OPTION_NUMBERING] = 0;

        std::vector<idx_t> part;

        run_metis_and_assign_comp(graph, component_array, part);
        for(int i=0; i<metis_setup.num_partition; i++)
        {
            std::string complex(COMPLEX_COMPONENT);
            component_type.push_back(complex + std::to_string(num_complex_comp));
            complex_comp_indicator.push_back(true);
        }

        // log component
        comps_num_edges.resize(curr_component_num,0);
        eip_t eip = boost::edges(graph);
        for(ei_t ei=eip.first; ei!=eip.second; ei++)
        {
            vd_t source = boost::source(*ei, graph);
            vd_t target = boost::target(*ei, graph);
            if(component_array[int(source)] == component_array[int(target)] && component_array[int(source)] != IMPOSSIBLE_COMP_NUM)
            {
                //shc_log_info(shc_logname, "component_array[int(source)] %d\n", component_array[int(source)]);
                //shc_log_info(shc_logname, "component_array[int(target)] %d\n", component_array[int(target)]);
                comps_num_edges[component_array[int(source)]] ++;
            }
        }

        if(metis_setup.is_multiple_partition)
        {
            shc_log_info(shc_logname, "Starting repartition\n");
            add_weight_to_cut_and_update_graph(graph, &part);
            create_metis_array_input(graph);
            run_metis_and_assign_comp(graph, component_array_aux, part);

            for(int i=0; i<metis_setup.num_partition; i++)
            {
                std::string re_partition(RE_PARTITION);
                std::string complex(COMPLEX_COMPONENT);
                component_type.push_back(re_partition+complex+std::to_string(num_complex_comp));
                complex_comp_indicator.push_back(true);
            }
            shc_log_info(shc_logname, "Finish repartition\n");

            // log component
            comps_num_edges.resize(curr_component_num,0);
            eip_t eip = boost::edges(graph);
            for(ei_t ei=eip.first; ei!=eip.second; ei++)
            {
                vd_t source = boost::source(*ei, graph);
                vd_t target = boost::target(*ei, graph);
                if(component_array_aux[int(source)] == component_array_aux[int(target)] && component_array_aux[int(source)] != IMPOSSIBLE_COMP_NUM)
                {
                    //shc_log_info(shc_logname, "component_array[int(source)] %d\n", component_array_aux[int(source)]);
                    //shc_log_info(shc_logname, "component_array[int(target)] %d\n", component_array[int(target)]);
                    comps_num_edges[component_array_aux[int(source)]] ++;
                }
            }
        }
        num_complex_comp++;


        //std::string contig_file = setting.local_files.output_path + "/contig_file";
        //log_contigs(contig_file);
    }
#ifdef LOG_CONTIG_GRAPH
    shc_log_info(shc_logname, "end\n");
#endif

}

void Contig_graph_handler::run_metis_and_assign_comp(graph_t & graph, std::vector<comp_num_t> & array,
                                                std::vector<idx_t> & part)
{
    idx_t num_vertices = boost::num_vertices(graph);
    part.resize(num_vertices*2);
    /*
    shc_log_info(shc_logname, "hello0\n");
    shc_log_info(shc_logname, "num of edges %d\n", boost::num_edges(graph));
    shc_log_info(shc_logname, "metis_setup.num_vertices %d\n", metis_setup.num_vertices);
    shc_log_info(shc_logname, "metis_setup.ncon %d\n", metis_setup.ncon);
    shc_log_info(shc_logname, "metis_input.xadj size %d\n", metis_input.xadj.size());
    shc_log_info(shc_logname, "metis_input.adjncy. size %d\n", metis_input.adjncy.size());
    shc_log_info(shc_logname, "metis_input.adjwgt. size %d\n", metis_input.adjwgt.size());
    shc_log_info(shc_logname, "metis_setup.num_partition  %d\n", metis_setup.num_partition);

    for(int i=0; i<metis_input.xadj.size(); i++)
        shc_log_info(shc_logname, "xadj %d\n", metis_input.xadj[i]);
    for(int i=0; i<metis_input.adjncy.size(); i++)
        shc_log_info(shc_logname, "adjncy %d, adjwgt %d\n", metis_input.adjncy[i],
                                                    metis_input.adjwgt[i]);
    shc_log_info(shc_logname, "part. size %d\n", part.size());
    */
    //shc_log_info(shc_logname, "before running metis\n");
    int ret = METIS_PartGraphKway(
          &metis_setup.num_vertices, &metis_setup.ncon,
          &metis_input.xadj.at(0), &metis_input.adjncy.at(0),
          NULL, NULL, &metis_input.adjwgt.at(0), &metis_setup.num_partition,
          NULL, NULL, metis_setup.options, &metis_setup.objval, &part.at(0));
    //shc_log_info(shc_logname, "after running metis\n");

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
        //shc_log_info(shc_logname, "contig index %u has comp %u\n",
        //        graph[static_cast<vd_t>(i)].contig_index,
        //        curr_component_num + part.at(i));
        array[graph[static_cast<vd_t>(i)].contig_index] =
        curr_component_num + part.at(i);
    }
    //log_comp_content();

    curr_component_num += metis_setup.num_partition;
}

void Contig_graph_handler::add_weight_to_cut_and_update_graph(graph_t & graph, std::vector<idx_t> * part)
{
#ifdef LOG_CONTIG_GRAPH
    shc_log_info(shc_logname, "Start readjust graph\n");
#endif
    typename boost::graph_traits < graph_t >::adjacency_iterator vi, vi_end;
    std::pair<ed_t, bool> aer;
    std::pair<vi_t, vi_t> vip = boost::vertices(graph);
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
#ifdef LOG_CONTIG_GRAPH
    shc_log_info(shc_logname, "Finish readjust graph\n");
#endif
}

void Contig_graph_handler::
get_connected_contig_index_without_graph( Kmer_counter_map_iterator & it,
      std::vector<contig_num_t> & contig_stack,
      char left_letter, char right_letter)
{
    contig_num_t contigA_index, contigB_index;
    uint64_t next_byte = 0;
    if (it == kh->kmer_counter.end())
    {
        shc_log_error("contig contains unknown kmer\n");
        exit(1);
    }
    contigA_index = it->second.contig;
    // for that kmer, check what extension available to it

    for(uint8_t i=0; i< BIT_PER_BYTE; i++)
    {
        total_bit_set_query++;
        if(IS_BIT_I_SET(it->second.info, i))
        {

            if(i<PREFIX_OFFSET)
            {
                if(left_letter == num_to_char[i]) // already in the contig
                    continue;
                next_byte = prepend_byte(&(it->first), i, kmer_length);
            }
            else
            {
                if(right_letter == num_to_char[i-PREFIX_OFFSET])
                    continue;
                next_byte = append_byte(&(it->first), i-PREFIX_OFFSET, kmer_length);
            }
            total_query++;
            Kmer_counter_map_iterator kmer_it = kh->kmer_counter.find(next_byte);
            if (kmer_it != kh->kmer_counter.end())
            {
                total_update++;
                contigB_index = kmer_it->second.contig;
                if (contigA_index != contigB_index)
                {

                    conn_contig_set.set_conn_contig(contigB_index);

                    if(explorable_contig_set.is_explorable(contigB_index))
                    {
                        //if contigB_index is found
                        contig_stack.push_back(contigB_index);
                        explorable_contig_set.set_explored(contigB_index);
                    }
                }
            }
        }
    }
}


void Contig_graph_handler::
get_connected_contig_index( Kmer_counter_map_iterator & it, graph_t & graph,
     Contig_vertex_map & contig_vertex_map, std::vector<contig_num_t> & contig_stack,
     vd_t curr_vd, char left_letter, char right_letter)
{
    contig_num_t contigA_index, contigB_index;
    uint64_t next_byte = 0;
    if (it == kh->kmer_counter.end())
    {
        shc_log_error("contig contains unknown kmer\n");
        exit(1);
    }
    contigA_index = it->second.contig;
    // for that kmer, check what extension available to it

    for(uint8_t i=0; i< BIT_PER_BYTE; i++)
    {
        total_bit_set_query++;
        if(IS_BIT_I_SET(it->second.info, i))
        {

            if(i<PREFIX_OFFSET)
            {
                if(left_letter == num_to_char[i]) // already in the contig
                    continue;
                next_byte = prepend_byte(&(it->first), i, kmer_length);
            }
            else
            {
                if(right_letter == num_to_char[i-PREFIX_OFFSET])
                    continue;
                next_byte = append_byte(&(it->first), i-PREFIX_OFFSET, kmer_length);
            }
            //decode_kmer(base, &(it->first), kmer_length);
            //decode_kmer(next_base, &next_byte, kmer_length);
            //shc_log_info(shc_logname, "%s has %u bit set, and next kmer %s \n",
            //                                            base, i, next_base);

            //in case the matched kmer is deleted since that does not belong
            //to any contig or that contig is deleted
            total_query++;
            Kmer_counter_map_iterator kmer_it = kh->kmer_counter.find(next_byte);
            if (kmer_it != kh->kmer_counter.end())
            {
                total_update++;
                contigB_index = kmer_it->second.contig;
                if (contigA_index != contigB_index)
                {
                    //shc_log_info(shc_logname, "Find new contig %u, before increment edge with me %u\n", contigB_index, contigA_index);
                    //shc_log_info(shc_logname, "%u, at %d\n", curr_contig, conn_point);
                    update_graph(contig_vertex_map, graph, curr_vd, contigB_index);

                    if(explorable_contig_set.is_explorable(contigB_index))
                    {
                        //if contigB_index is found
                        contig_stack.push_back(contigB_index);
                        explorable_contig_set.set_explored(contigB_index);
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
void Contig_graph_handler::
update_graph(Contig_vertex_map & contig_vertex_map, graph_t & graph, vd_t vd_i, contig_num_t j)
{
    vd_t vd_j, it_j;
    ;
    //check if we need to add j
    if(!contig_vertex_map.get_contig_vd(j, vd_j))
    {
        //shc_log_info(shc_logname, "without %u\n", j);
        vd_j = boost::add_vertex(graph);
        contig_vertex_map.set_contig_vd(j, vd_j) ;
        graph[vd_j] = bundled_contig_index(j);
#ifdef LOG_CONTIG_GRAPH
        shc_log_info(shc_logname, "Added contig %u with vd %u\n", j, vd);
#endif
    }

    // note: in boost graph if it is an undirectedS, and edge between vd_i and
    // vd_j exists, then check on vd_j and vd_i shows true
    std::pair<ed_t, bool> aer = boost::edge(vd_i, vd_j, graph);
    if(aer.second)
    {
        (graph[aer.first].weight)++;

#ifdef LOG_CONTIG_GRAPH
        //shc_log_info(shc_logname, "edge between %u %u : weight increased to %u\n",
        //        i, j ,graph[aer.first].weight);
#endif
    }
    else //add edge
    {
        aer = boost::add_edge(vd_i, vd_j, graph);
        graph[aer.first] = bundled_weight(1);
#ifdef LOG_CONTIG_GRAPH
        //shc_log_info(shc_logname, "Edge between %u %u\n", i, j);
#endif
    }
}

void Contig_graph_handler::create_metis_array_input(graph_t & graph)
{
#ifdef LOG_CONTIG_GRAPH
    shc_log_info(shc_logname, "Start creating metis array input for component\n",
                                                curr_component_num);
#endif
    //shc_log_info(shc_logname, "Start creating metis array input for component %d\n",
    //                                            curr_component_num);
    //shc_log_info(shc_logname, "num_edges %d, num nodes %d\n",
    //                    boost::num_edges(graph), boost::num_vertices(graph));
    metis_input.metis_reset_data(boost::num_edges(graph), boost::num_vertices(graph));

    typedef boost::graph_traits< graph_t >::out_edge_iterator out_ei_t;
    typedef std::pair<out_ei_t, out_ei_t> out_eip_t;
    contig_edge_weight_t edge_weight;

    std::pair<vi_t, vi_t> vip = boost::vertices(graph);
    int i=0, j=1;
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        vd_t vd_source = *it;
        out_eip_t out_eip = boost::out_edges(vd_source,graph);
        for(out_ei_t out_ei=out_eip.first; out_ei!=out_eip.second; out_ei++)
        {
            vd_t vd_target = boost::target(*out_ei, graph);
            edge_weight = (graph[*out_ei].weight)/2;
            metis_input.adjncy[i] = vd_target;
            metis_input.adjwgt[i] = edge_weight;
            i++;
        }
        metis_input.xadj[j++] = i;
    }
    //shc_log_info(shc_logname, "i: %d, 2edge: %d\n", i, 2*boost::num_edges(graph));
    //shc_log_info(shc_logname, "j: %d, node: %d\n", j, boost::num_vertices(graph));
    assert(i==2*boost::num_edges(graph));
    assert(j==boost::num_vertices(graph)+1);
    shc_log_info(shc_logname, "Finish\n");
    //log_metis_input_data();
#ifdef LOG_CONTIG_GRAPH
    shc_log_info(shc_logname, "Finish creating metis array input for component, %d %d %d\n",
                        metis_input.xadj.size(), metis_input.adjncy.size(), metis_input.adjwgt.size());
#endif
}
/*
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

    std::pair<vi_t, vi_t> vip = boost::vertices(graph);

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
            edge_weight = graph[aer.first].weight;
            outfile << vd << " " << edge_weight << " ";    //without offset
        }
        outfile << std::endl;
    }

    //log_metis_input_data();
    outfile.close();
    shc_log_info(shc_logname, "Finish making METIS format file\n");
}
*/
/*
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
*/

void Contig_graph_handler::dump_component_array(std::string & filename)
{
    shc_log_info(shc_logname, "Start dumping: %u component \n", curr_component_num);
    std::ofstream outfile(filename);
    outfile << curr_component_num
      << " components, this number also denotes unmmapped contig to components, each row is a contig"
      << std::endl;

    //int comp_type_i = 0;
    //for(std::vector<std::string>::iterator it=component_type.begin();
    //                                       it!=component_type.end(); it++)
    //{
    //    shc_log_info(shc_logname, "%d: comp type %s\n",(comp_type_i++) ,(*it).c_str());
    //}

    //int i=0;
    for(int i=0; i<component_array.size(); i++)
    {
        if(!metis_setup.is_multiple_partition)
        {
            outfile << component_array[i] << "\t" << (component_type[component_array[i]]) << std::endl;
        }
        else
        {
            if(component_array_aux[i] == IMPOSSIBLE_COMP_NUM)
                component_array_aux[i] = curr_component_num;

            outfile << (comp_num_t)component_array[i] << "\t"
                    << (comp_num_t)component_array_aux[i] << "\t"
                    << (component_type[component_array[i]]) << std::endl;
        }
        //shc_log_info(shc_logname, "%d   %u\n", i++, *it);
    }
    outfile.close();
    shc_log_info(shc_logname, "Finish log_component_array\n");

    std::ofstream comp_type_writer(lf.output_path + "/comp_types.log");
    comp_type_writer << curr_component_num
      << " components, c for complex, s for simple"
      << std::endl;
    for (int i=0; i<curr_component_num; i++)
    {
        comp_type_writer << component_type[i] << std::endl;
    }
    comp_type_writer.close();
}


void Contig_graph_handler::
update_comp_count(std::vector<comp_num_t> & components,
                  std::vector<comp_num_t> & counts, bool is_forward,
                  std::vector<comp_num_t> & comp_array,  char * seq,  int len)
{
    int interval = (len-kmer_length+1)/num_feature;
    uint64_t byte;
    for(int i=0; i<num_feature; i++)
    {
        if(is_forward)
            encode_kmer(seq+i*interval, &byte, kmer_length);
        else
        {
            encode_reverse_kmer(seq+len-i*interval-1, &byte, kmer_length);
            complement_num(&byte, kmer_length);
        }


        Kmer_counter_map_iterator it = kh->kmer_counter.find(byte);
        if(it != kh->kmer_counter.end())
        {
            contig_num_t contig_num_for_kmer = it->second.contig;
            comp_num_t comp_num_for_contig = comp_array[contig_num_for_kmer];
            if(comp_num_for_contig < curr_component_num) // if equal, then it
            {
                std::vector<comp_num_t>::iterator it = std::find(
                      components.begin(), components.end(), comp_num_for_contig);
                if(it==components.end())
                {
                    components.push_back(comp_num_for_contig);
                    counts.push_back(1);
                }
                else
                {
                    counts[it-components.begin()]++;
                }
            }
        }
    }
}

void Contig_graph_handler::
update_comp_map(std::map<comp_num_t, int> & comp_map, bool is_forward,
                std::vector<comp_num_t> & comp_array,  char * seq,  int len,
                kmer_count_t * kmer_counts, int & num_kmer)
{
    int interval = (len-kmer_length+1)/num_feature;
    num_kmer = 0;
    //int num_featureed = 0;
    uint64_t byte;

    //bool is_set_min_kmer = false;
    //kmer_count_t min_kmer_count;

    //std::string seq_str(seq, len);
    //shc_log_info(shc_logname, "seq %s\n", seq_str.c_str());
    //shc_log_info(shc_logname, "is_forward %s\n", (is_forward)?("yes"):("no"));
    //shc_log_info(shc_logname, "num_feature %s\n", interval);

    //char de_kmer[33];
    //de_kmer[kmer_length] = '\0';

    //while(offset < interval)
    //{
        for(int i=0; i<num_feature; i++)
        {
            if(is_forward)
            {
                int forward_i = i*interval;

                //shc_log_info(shc_logname, "forward check\n");
                if(!encode_kmer_check_base(seq+forward_i, &byte, kmer_length))
                {
                    comp_map.clear();
                    return;
                    //continue;
                }
                //shc_log_info(shc_logname, "forward check F\n");
            }
            else
            {
                //shc_log_info(shc_logname, "reverse check\n");
                int rc_i =  len-i*interval-kmer_length-1; //i*interval

                if(!encode_reverse_kmer_check_base(seq+rc_i, &byte, kmer_length))
                {
                    comp_map.clear();
                    return;
                }
                //continue;

                complement_num(&byte, kmer_length);
                //shc_log_info(shc_logname, "reverse check F\n");
            }

            //decode_kmer(de_kmer, &byte, kmer_length);
            //shc_log_info(shc_logname, "kmer %s\n", de_kmer);

            Kmer_counter_map_iterator it = kh->kmer_counter.find(byte);
            if(it != kh->kmer_counter.end())
            {
                contig_num_t contig_num_for_kmer = it->second.contig;
                comp_num_t comp_num_for_contig = comp_array[contig_num_for_kmer];
                if(comp_num_for_contig < curr_component_num) // if equal, then it
                {
                    comp_map[comp_num_for_contig]++;
                }
                //kmer_count_t kmer_count = ;
                //if(!is_set_min_kmer)
                //{
                //    is_set_min_kmer = true;
                //    min_kmer_count = kmer_count;
                //}
                //else
                //{
                //    if (kmer_count < min_kmer_count )
                //        min_kmer_count = kmer_count;
                //}
                assert(NUM_FEAT_RESERVED > num_kmer);
                kmer_counts[num_kmer++] = get_count(it->second.count, compress_ratio);

                //kmer_counts.push_back(get_count(it->second.count, compress_ratio));

                //avg_count += get_count(it->second.count, compress_ratio);
                //shc_log_info(shc_logname, "count %u\n", get_count(it->second.count, compress_ratio));
                //num_valid_test ++;
            }
            else
            {
                //num_kmer = 0;
                comp_map.clear();
                return;
            }
        }
    //    offset ++;
    //    if (num_valid_test >= num_test)
    //        break;
    //}
    // it is ok to miss, since they are erroreous kmer, the ream
    /*
    if (kmer_counts.size() > 0)
    {
        if(num_feature%2 == 1)
        {
            int median_i = num_feature/2;
            const auto median_it  = kmer_counts.begin()+median_i;
            std::nth_element (kmer_counts.begin(),
                            median_it, kmer_counts.end());
            avg_count = *median_it;
        }
        else
        {
            const auto median_it1 = kmer_counts.begin() + kmer_counts.size() / 2 - 1;
            const auto median_it2 = kmer_counts.begin() + kmer_counts.size() / 2;

            std::nth_element(kmer_counts.begin(), median_it1 , kmer_counts.end());
            const auto e1 = *median_it1;

            std::nth_element(kmer_counts.begin(), median_it2 , kmer_counts.end());
            const auto e2 = *median_it2;

            avg_count = (e1 + e2) / 2;
        }
    }
    */
    //if(is_set_min_kmer)
    //    avg_count = min_kmer_count;
    //else
    //    avg_count = 0;
    //shc_log_info(shc_logname, "avg_count %f\n", avg_count);
}

void Contig_graph_handler::
update_comp_set(std::set<comp_num_t> & comp_set, bool is_forward,
        std::vector<comp_num_t> & comp_array, char * seq,  int len)
{
    int interval = (len-kmer_length+1)/num_feature;
    uint64_t byte;
    for(int i=0; i<num_feature; i++)
    {
        if(is_forward)
            encode_kmer(seq+i*interval, &byte, kmer_length);
        else
        {
            encode_reverse_kmer(seq+len-i*interval-1, &byte, kmer_length);
            complement_num(&byte, kmer_length);
        }


        Kmer_counter_map_iterator it = kh->kmer_counter.find(byte);
        if(it != kh->kmer_counter.end())
        {
            contig_num_t contig_num_for_kmer = it->second.contig;
            comp_num_t comp_num_for_contig = comp_array[contig_num_for_kmer];
            if(comp_num_for_contig < curr_component_num) // if equal, then it
            {
                comp_set.insert(comp_num_for_contig);
            }

        }
    }
}

bool Contig_graph_handler::
decide_best_comp(std::map<comp_num_t, int> & comp_count, comp_num_t & best_comp)
{
    if(!comp_count.empty())
    {
        auto max_iter = std::max_element(comp_count.begin(), comp_count.end(),
                [](const std::pair<comp_num_t, int> & p1,
                   const std::pair<comp_num_t, int> & p2 )
                {return p1.second < p2.second;});

        best_comp =  max_iter->first;
        return true;
    }
    best_comp = curr_component_num;
    return false;
}

void Contig_graph_handler::assign_paired_read_to_components_mmap(
                        std::string & read_1_path, std::string & read_2_path)
{
    shc_log_info(shc_logname, "Start assigning paired read to component\n");
    //creating dir
    typedef boost::filesystem::path path_t;
    std::string dir_path = lf.output_components_read_dir;

    std::string file_prefix = dir_path + lf.comp_read_prefix;
    std::string file_p1_suffix("_p1");
    std::string file_p2_suffix("_p2");

    struct Single_dumper & dumper_p1 = fasta_dumper.pair1_dumper;
    struct Single_dumper & dumper_p2 = fasta_dumper.pair2_dumper;

    Block_timer p_timer;
    start_timer(&p_timer);
    dumper_p1.setup_all_mmap_dst_files(file_prefix, file_p1_suffix ,
                curr_component_num);
    dumper_p2.setup_all_mmap_dst_files(file_prefix, file_p2_suffix ,
                curr_component_num);

    //std::cout << "dst file setup time takes " << std::endl;
    //stop_timer(&p_timer);
    shc_log_info(shc_logname, "dst file setup time takes \n");
    log_stop_timer(&p_timer);

    start_timer(&p_timer);
    dumper_p1.setup_init_mmap_src_file(read_1_path); //mmap_src_file(read_1_path);
    dumper_p2.setup_init_mmap_src_file(read_2_path); //mmap_src_file
    //std::cout << "src file setup time takes " << std::endl;
    //stop_timer(&p_timer);
    shc_log_info(shc_logname, "src file setup time takes \n");
    log_stop_timer(&p_timer);

    uint64_t num_line = 0;
    int64_t num_read =0 ;
    int64_t last_num = 0;

    char * header_ptr_p1;
    char * header_ptr_p2;
    char * seq_ptr_p1;
    char * seq_ptr_p2;
    int seq_len_p1, seq_len_p2;

    Block_timer read_prog_timer;
    Block_timer all_time;
    start_timer(&read_prog_timer);
    start_timer(&all_time);

    {
        std::string message = "Assign pair-reads to components\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t progress_step = num_pair_read_file/100;

        while(dumper_p1.get_src_line(header_ptr_p1, seq_ptr_p1, seq_len_p1) &&
              dumper_p2.get_src_line(header_ptr_p2, seq_ptr_p2, seq_len_p2))
        {
            assign_pair_read_to_file_mmap(dumper_p1, dumper_p2,
                                          header_ptr_p1, header_ptr_p2,
                                          seq_ptr_p1,    seq_ptr_p2,
                                          seq_len_p1,    seq_len_p2);
            num_read++;
            if((num_read -last_num) > progress_step)
            {
                int percentage = (100 * num_read)/num_pair_read_file;
                progress.write(static_cast<double>(num_read) /
                               static_cast<double>(num_pair_read_file));

                shc_log_info(shc_logname, "[ %u%] %u read out of %u is processed\n",
                        percentage, num_read, (num_pair_read_file));
                log_stop_timer(&read_prog_timer);
                start_timer(&read_prog_timer);
                last_num = num_read;
            }
        }
    }

    shc_log_info(shc_logname, "finish dump pair read\n");
    //std::cout << "finish dump pair read" << std::endl;
    //stop_timer(&all_time);
    log_stop_timer(&all_time);
}

void Contig_graph_handler::assign_paired_read_to_components(
                    std::string & read_1_path, std::string & read_2_path)
{
    shc_log_info(shc_logname, "Start assigning paired read to component\n");
    //creating dir
    typedef boost::filesystem::path path_t;
    std::string dir_path = lf.output_components_read_dir;
    path_t comp_read_path(dir_path);
    if(boost::filesystem::exists(comp_read_path))
        empty_directory(comp_read_path, lf.output_path);
    else
        add_directory(comp_read_path);

    //create and open files
    std::vector<std::shared_ptr<std::ofstream> > files1;
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        std::string file_path = dir_path +
                        (lf.comp_read_prefix + std::to_string(i)+ "_p1");
        std::shared_ptr<std::ofstream> file(new std::ofstream);
        file->open(file_path.c_str());
        files1.push_back(file);
    }
    std::vector<std::shared_ptr<std::ofstream> > files2;
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        std::string file_path = dir_path +
                        (lf.comp_read_prefix + std::to_string(i) + "_p2");
        std::shared_ptr<std::ofstream> file(new std::ofstream);
        file->open(file_path.c_str());
        files2.push_back(file);
    }

    shc_log_info(shc_logname, "Created contig files\n");
    std::ifstream file1_reader(read_1_path);
    std::ifstream file2_reader(read_2_path);
    std::string line1, header1, sequence1;
    std::string line2, header2, sequence2;

    while (std::getline(file1_reader, line1) && std::getline(file2_reader, line2))
    {
        if(line1[0] == '>' && line2[0] == '>')
        {
            if(!sequence1.empty() && !sequence2.empty())
            {
                assign_pair_read_to_file(header1, sequence1,
                                         header2, sequence2,
                                         files1, files2);
                // here we consider there two reads in case of double strand
                // read is only assigned to its best comp
                if(setting.is_double_stranded)
                {
                    std::reverse(sequence1.begin(), sequence1.end());
                    std::reverse(sequence2.begin(), sequence2.end());
                    // note we swap files name, since it is direction is changed
                    assign_pair_read_to_file(header1, sequence1,
                                         header2, sequence2,
                                         files2, files1);
                }
                sequence1.clear(), sequence2.clear();
            }
            header1 = line1;
            header2 = line2;
        }
        else //line is sequence
        {
            sequence1 += line1;
            sequence2 += line2;
        }
    }

    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        (*files1.at(i)).close();
        (*files2.at(i)).close();
    }

    std::vector<std::shared_ptr<std::ofstream> >().swap(files1);
    std::vector<std::shared_ptr<std::ofstream> >().swap(files2);
    shc_log_info(shc_logname, "Finish assigning read to component\n");
}

void Contig_graph_handler::
assign_pair_read_to_file(std::string & header1, std::string & sequence1,
                  std::string & header2, std::string & sequence2,
                  std::vector<std::shared_ptr<std::ofstream> > & files1,
                  std::vector<std::shared_ptr<std::ofstream> > & files2)
{
#ifdef LOG_CONTIG_GRAPH
    shc_log_info(shc_logname, "Start assigning a pair of read\n");
#endif

    std::set<comp_num_t> comp_set;

    update_comp_set(comp_set, true, component_array, &sequence1.at(0), sequence1.size());
    update_comp_set(comp_set, true, component_array, &sequence2.at(0), sequence2.size());
    if (metis_setup.is_multiple_partition)
    {
        update_comp_set(comp_set, true, component_array_aux, &sequence1.at(0), sequence1.size());
        update_comp_set(comp_set, true, component_array_aux, &sequence2.at(0), sequence2.size());
    }
    for(std::set<comp_num_t>::iterator it =comp_set.begin();
                                       it!=comp_set.end(); it++)
    {
        *(files1.at(*it)) << header1 << std::endl;
        *(files1.at(*it)) << sequence1 << std::endl;
        *(files2.at(*it)) << header2 << std::endl;
        *(files2.at(*it)) << sequence2 << std::endl;
    }
#ifdef LOG_CONTIG_GRAPH
    shc_log_info(shc_logname, "Finish assigning a pair of read\n");
#endif
}


/*  // choose all
void Contig_graph_handler::
assign_single_read_to_file_mmap(char * header_ptr, char * seq_ptr, int seq_len)
{
    //shc_log_info(shc_logname, "start\n");
    //std::vector<comp_num_t> components;
    //std::vector<comp_num_t> counts;
    std::map<comp_num_t, int> comp_count_map;
    update_comp_map(comp_count_map, true, component_array, seq_ptr, seq_len);
    //update_comp_count(components, counts, true, component_array, seq_ptr, seq_len);
    if (metis_setup.is_multiple_partition)
        update_comp_map(comp_count_map, true, component_array_aux, seq_ptr, seq_len);
    for(std::map<comp_num_t, int>::iterator it=comp_count_map.begin();
                                          it!= comp_count_map.end(); it++)
    {
        comp_num_t comp_i = it->first;
        uint64_t len = seq_ptr-header_ptr+seq_len+1;
        if(dumper.is_increase_file_sz(comp_i, len))
        {
            mmap_files.remap_file_increment_size(comp_i);
        }
        mmap_files.mmap_write(comp_i, len, header_ptr);
    }

    //shc_log_info(shc_logname, "before double stranded\n");
    if(is_double_stranded)
    {
        std::map<comp_num_t, int> ds_comp_map;
        update_comp_map(ds_comp_map, false, component_array, seq_ptr, seq_len);
        if (metis_setup.is_multiple_partition)
            update_comp_map(ds_comp_map, false, component_array_aux, seq_ptr, seq_len);
        for(std::map<comp_num_t, int>::iterator it=ds_comp_map.begin();
                                              it!= ds_comp_map.end(); it++)
        {
            comp_num_t comp_i = it->first;
            uint64_t len = seq_ptr-header_ptr+seq_len+1;
            if(mmap_files.is_increase_file_sz(comp_i, len))
            {
                mmap_files.remap_file_increment_size(comp_i);
            }
            mmap_files.mmap_write(comp_i, len, header_ptr);
        }
    }
    //shc_log_info(shc_logname, "end\n");
}
*/

void Contig_graph_handler::
assign_pair_read_to_file_mmap_helper(struct Single_dumper & mfs_p1,
             struct Single_dumper & mfs_p2, char * header_ptr_p1,
             char * header_ptr_p2, char * seq_ptr_p1, char * seq_ptr_p2,
             int seq_len_p1, int seq_len_p2, bool is_forward)
{
    //if(is_assign_best)
    //{
        int num_kmer1, num_kmer2;
        std::map<comp_num_t, int> comp_count_map_p1;
        std::map<comp_num_t, int> comp_count_map_p2;
        // read p1 update and dump
        update_comp_map(comp_count_map_p1, is_forward, component_array,
                        seq_ptr_p1, seq_len_p1, (kmer_count_t*)pair_1_kmer_counts, num_kmer1);
        update_comp_map(comp_count_map_p2, is_forward, component_array,
                        seq_ptr_p2, seq_len_p2, (kmer_count_t*)pair_2_kmer_counts, num_kmer2);
        comp_num_t best_comp_p1, best_comp_p2;

        if(decide_best_comp(comp_count_map_p1, best_comp_p1))
        {
            assign_pair_best_comp(mfs_p1, mfs_p2,
               best_comp_p1, seq_ptr_p1, seq_ptr_p2, header_ptr_p1,
               header_ptr_p2, seq_len_p1, seq_len_p2, is_forward, num_kmer1, num_kmer2);
        }

        if(decide_best_comp(comp_count_map_p2, best_comp_p2))
        {
            if(best_comp_p1 != best_comp_p2)
            {
                assign_pair_best_comp(mfs_p1, mfs_p2,
                   best_comp_p2, seq_ptr_p1, seq_ptr_p2, header_ptr_p1,
                   header_ptr_p2, seq_len_p1, seq_len_p2, is_forward, num_kmer1, num_kmer2);
            }
        }

        if (metis_setup.is_multiple_partition)
        {

            int re_num_kmer1, re_num_kmer2;
            comp_count_map_p1.clear();
            update_comp_map(comp_count_map_p1, is_forward,
                          component_array_aux, seq_ptr_p1, seq_len_p1,
                          (kmer_count_t*)pair_1_kmer_counts, re_num_kmer1);

            comp_count_map_p2.clear();
            update_comp_map(comp_count_map_p2, is_forward,
                        component_array_aux, seq_ptr_p2, seq_len_p2,
                        (kmer_count_t*)pair_2_kmer_counts, re_num_kmer2);

             // is complex component
            comp_num_t re_best_comp_p1 = curr_component_num;
            if((best_comp_p1!=curr_component_num && complex_comp_indicator[best_comp_p1]))
            {
                if(decide_best_comp(comp_count_map_p1, re_best_comp_p1))
                {
                    assign_pair_best_comp(mfs_p1, mfs_p2,
                         re_best_comp_p1, seq_ptr_p1, seq_ptr_p2, header_ptr_p1,
                         header_ptr_p2, seq_len_p1, seq_len_p2, is_forward,
                         re_num_kmer1, re_num_kmer2);
                }
            }

            if((best_comp_p2!=curr_component_num && complex_comp_indicator[best_comp_p2]))
            {
                comp_num_t re_best_comp_p2 = curr_component_num;
                if(decide_best_comp(comp_count_map_p2, re_best_comp_p2))
                {
                    if(re_best_comp_p1 != re_best_comp_p2)
                    {
                        assign_pair_best_comp(mfs_p1, mfs_p2,
                           re_best_comp_p2, seq_ptr_p1, seq_ptr_p2, header_ptr_p1,
                           header_ptr_p2, seq_len_p1, seq_len_p2, is_forward,
                           re_num_kmer1, re_num_kmer2);
                    }
                }
            }
        }
    //}
    //else
    //{
    //    int num_kmer1, num_kmer2;
    //    std::map<comp_num_t, int> comp_count_map;
    //    update_comp_map(comp_count_map, is_forward, component_array,
    //                    seq_ptr_p1, seq_len_p1, num_kmer1);
    //    update_comp_map(comp_count_map, is_forward, component_array,
    //                    seq_ptr_p2, seq_len_p2, num_kmer2);
    //    if (metis_setup.is_multiple_partition)
    //    {
    //        update_comp_map(comp_count_map, is_forward, component_array_aux,
    //                        seq_ptr_p1, seq_len_p1, num_kmer1);
    //        update_comp_map(comp_count_map, is_forward, component_array_aux,
    //                        seq_ptr_p2, seq_len_p2, num_kmer2);
    //    }

    //    assign_pair_all_comp(mfs_p1, mfs_p2, comp_count_map,
    //            seq_ptr_p1, seq_ptr_p2, header_ptr_p1, header_ptr_p2,
    //            seq_len_p1, seq_len_p2, is_forward);
    //}
}


comp_num_t Contig_graph_handler::
assign_pair_best_comp(struct Single_dumper & mfs_p1,
             struct Single_dumper & mfs_p2, comp_num_t comp_i,
                char * seq_ptr_p1, char * seq_ptr_p2, char * header_ptr_p1,
                char * header_ptr_p2, int seq_len_p1, int seq_len_p2, bool is_forward,
                int num_kmer1, int num_kmer2)
{
    if( read_sampler.is_store_all_reads_and_features ||
        read_sampler.decide_to_keep_read(comp_i, pair_1_kmer_counts, num_kmer1) ||
        read_sampler.decide_to_keep_read(comp_i, pair_2_kmer_counts, num_kmer2) )
    {
        if(is_forward)
        {
            mfs_p1.write_fasta(comp_i, seq_ptr_p1, header_ptr_p1, seq_len_p1);
            mfs_p2.write_fasta(comp_i, seq_ptr_p2, header_ptr_p2, seq_len_p2);
        }
        else
        {
            mfs_p1.write_reverse_fasta(comp_i, seq_ptr_p2, header_ptr_p2, seq_len_p2);
            mfs_p2.write_reverse_fasta(comp_i, seq_ptr_p1, header_ptr_p1, seq_len_p1);
        }
        if(read_sampler.is_store_all_reads_and_features)
        {
            //read_features_dumper.dump_read_prob(comp_i,
            //    std::min(std::max(prob_to_sample_1, prob_to_sample_2), 1.0) );
            paired_read_features_dumper.dump_paired_read_features(comp_i,
                        pair_1_kmer_counts, num_kmer1,
                        pair_2_kmer_counts, num_kmer2);
        }

        return comp_i;
    }
    else
        return curr_component_num;
}

void Contig_graph_handler::
assign_pair_all_comp(struct Single_dumper & mfs_p1,
             struct Single_dumper & mfs_p2, std::map<comp_num_t, int> & comp_count,
    char * seq_ptr_p1, char * seq_ptr_p2, char * header_ptr_p1,
                char * header_ptr_p2, int seq_len_p1, int seq_len_p2, bool is_forward)
{
    for(std::map<comp_num_t, int>::iterator it=comp_count.begin();
                                          it!= comp_count.end(); it++)
    {
        comp_num_t comp_i = it->first;
        if(is_forward)
        {
            mfs_p1.write_fasta(comp_i, seq_ptr_p1, header_ptr_p1,
                seq_len_p1);
            mfs_p2.write_fasta(comp_i, seq_ptr_p2, header_ptr_p2,
                seq_len_p2);
        }
        else
        {
            mfs_p1.write_reverse_fasta(comp_i, seq_ptr_p2, header_ptr_p2, seq_len_p2);
            mfs_p2.write_reverse_fasta(comp_i, seq_ptr_p1, header_ptr_p1, seq_len_p1);
        }
    }
}

void Contig_graph_handler::
assign_pair_read_to_file_mmap(struct Single_dumper & mfs_p1,
             struct Single_dumper & mfs_p2, char * header_ptr_p1,
             char * header_ptr_p2, char * seq_ptr_p1, char * seq_ptr_p2,
             int seq_len_p1, int seq_len_p2)
{
    assign_pair_read_to_file_mmap_helper(mfs_p1,
                 mfs_p2, header_ptr_p1,
                 header_ptr_p2, seq_ptr_p1, seq_ptr_p2,
                 seq_len_p1, seq_len_p2, true);

    //shc_log_info(shc_logname, "before double stranded\n");
    if(is_double_stranded)
    {
        assign_pair_read_to_file_mmap_helper(mfs_p1,
                     mfs_p2, header_ptr_p1,
                     header_ptr_p2, seq_ptr_p1, seq_ptr_p2,
                     seq_len_p1, seq_len_p2, false);
    }
}

 //    choose best only
void Contig_graph_handler::
assign_single_read_to_file_mmap(struct Single_dumper & mfs,
                        char * header_ptr, char * seq_ptr, int seq_len)
{
    //shc_log_info(shc_logname, "assign_single_read_to_file_mmap\n");
    assign_single_read_to_file_mmap_helper(mfs, header_ptr, seq_ptr, seq_len, true);
    //char test[100000];
    //memset(test, 0, 100000);
    //memcpy(test, seq_ptr, seq_len);
    //shc_log_info(shc_logname, "%s\n", test);
    if(is_double_stranded)
    {
        assign_single_read_to_file_mmap_helper(mfs, header_ptr, seq_ptr, seq_len, false);
    }
    //shc_log_info(shc_logname, "finish_single_read_to_file_mmap\n");
}

comp_num_t Contig_graph_handler::
assign_single_best_comp(struct Single_dumper & mfs,
    comp_num_t comp_i, char * seq_ptr, char * header_ptr,
                                int seq_len, bool is_forward, int num_kmer)
{
    if( read_sampler.is_store_all_reads_and_features  ||
        read_sampler.decide_to_keep_read(comp_i, single_kmer_counts, num_kmer))
    {
        if (is_forward)
            mfs.write_fasta(comp_i, seq_ptr, header_ptr, seq_len);
        else
            mfs.write_reverse_fasta(comp_i, seq_ptr, header_ptr, seq_len);

        if(read_sampler.is_store_all_reads_and_features)
            single_read_features_dumper.dump_read_features(comp_i, single_kmer_counts, num_kmer);
        //char read_out[1000];
        //uint64_t len = seq_ptr-header_ptr;
        //memcpy(read_out, header_ptr, len);
        //read_out[len] = '\0';
        //shc_log_info(shc_logname, "%s", read_out);
        //memcpy(read_out, seq_ptr, seq_len);


        return comp_i;
    }
    else
        return curr_component_num;
}

void Contig_graph_handler::
assign_single_all_comp(struct Single_dumper & mfs,
    std::map<comp_num_t, int> & comp_count, char * seq_ptr, char * header_ptr,
                                        int seq_len, bool is_forward)
{
    for(std::map<comp_num_t, int>::iterator it=comp_count.begin();
                                          it!= comp_count.end(); it++)
    {
        comp_num_t comp_i = it->first;
        if (is_forward)
            mfs.write_fasta(comp_i, seq_ptr, header_ptr, seq_len);
        else
            mfs.write_reverse_fasta(comp_i, seq_ptr, header_ptr, seq_len);
    }
}

void Contig_graph_handler::
assign_single_read_to_file_mmap_helper(struct Single_dumper & mfs,
                        char * header_ptr, char * seq_ptr, int seq_len, bool is_forward)
{
    //if(is_assign_best)
    //{
        int num_tested_kmer;
        std::map<comp_num_t, int> comp_count_map;
        comp_num_t best_comp = curr_component_num;
        update_comp_map(comp_count_map, is_forward, component_array,
                        seq_ptr, seq_len, (kmer_count_t*)single_kmer_counts, num_tested_kmer);
        if(decide_best_comp(comp_count_map, best_comp))
        {
            assign_single_best_comp(mfs, best_comp, seq_ptr, header_ptr,
                                                    seq_len, is_forward, num_tested_kmer);
        }

        if(metis_setup.is_multiple_partition)
        {
            if(best_comp!=curr_component_num &&complex_comp_indicator[best_comp]) // is complex component
            {
                int re_num_tested_kmer;
                comp_count_map.clear();
                update_comp_map(comp_count_map, is_forward,
                        component_array_aux, seq_ptr, seq_len,
                        (kmer_count_t*)single_kmer_counts, re_num_tested_kmer);
                if(decide_best_comp(comp_count_map, best_comp))
                {
                    assign_single_best_comp(mfs, best_comp, seq_ptr,
                                     header_ptr, seq_len, is_forward, re_num_tested_kmer);
                }
            }
        }
    //}
    //else
    //{
    //    double avg_count;
    //    std::map<comp_num_t, int> comp_count_map;
    //    update_comp_map(comp_count_map, is_forward, component_array,
    //                    seq_ptr, seq_len, avg_count);
    //    if (metis_setup.is_multiple_partition)
    //    {
    //        update_comp_map(comp_count_map, is_forward, component_array_aux,
    //                        seq_ptr, seq_len, avg_count);
    //    }
    //    assign_single_all_comp(mfs, comp_count_map, seq_ptr, header_ptr,
    //                                        seq_len, is_forward);
    //}
}


void Contig_graph_handler::
assign_single_read_to_file(std::string & header, std::string & sequence,
             std::vector<std::shared_ptr<std::ofstream> > & files)
{
    std::set<comp_num_t> comp_set;

    update_comp_set(comp_set, true, component_array, &sequence.at(0), sequence.size());
    if (metis_setup.is_multiple_partition)
        update_comp_set(comp_set, true, component_array_aux, &sequence.at(0), sequence.size());
    for(std::set<comp_num_t>::iterator it=comp_set.begin(); it!= comp_set.end(); it++)
    {
        *(files.at(*it)) << header << std::endl;
        *(files.at(*it)) << sequence << std::endl;
    }
}

void Contig_graph_handler::
assign_reads_to_components(std::string input_read_path)
{
    start_timer(&cgh_timer);
    shc_log_info(shc_logname, "Start assigning read to component\n");
    //creating dir
    typedef boost::filesystem::path path_t;
    std::string dir_path = lf.output_components_read_dir;

    //create files
    std::vector<std::shared_ptr<std::ofstream> > files;
    for(comp_num_t i=0; i<curr_component_num; i++)
    {
        std::string file_path = dir_path +
                        (lf.comp_read_prefix + std::to_string(i));
        std::shared_ptr<std::ofstream> file(new std::ofstream);
        file->open(file_path.c_str());
        files.push_back(file);
    }
    shc_log_info(shc_logname, "Created contig files\n");

    //read and process files
    std::ifstream file_reader(input_read_path);
    std::string line;
    std::string header;
    std::string sequence;

    while (std::getline(file_reader, line))
    {
        if(line[0] == '>')
        {
            if(!sequence.empty())
            {
                assign_single_read_to_file(header, sequence, files);
                // here we consider there two reads in case of double strand,
                // read is only assigned to its best comp
                if(setting.is_double_stranded)
                {
                    std::reverse(sequence.begin(), sequence.end()); // in place
                    assign_single_read_to_file(header, sequence, files);
                }
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
    //std::cout << "finish dump reads " ;
    //stop_timer(&cgh_timer);
}

void Contig_graph_handler::
assign_reads_to_components_mmap(std::string input_read_path)
{
    start_timer(&cgh_timer);
    shc_log_info(shc_logname, "Start assigning read to component\n");
    //creating dir
    typedef boost::filesystem::path path_t;
    std::string dir_path = lf.output_components_read_dir;

    std::string file_prefix = dir_path + lf.comp_read_prefix;
    std::string file_suffix(""); //no suffix

    struct Single_dumper & dumper = fasta_dumper.single_dumper;
    dumper.setup_all_mmap_dst_files(file_prefix, file_suffix,
                curr_component_num);
    //dumper.mmap_src_file(input_read_path);

    dumper.setup_init_mmap_src_file(input_read_path);

    uint64_t last_num = 0;

    //printf("statbuf.st_size %ld\n", dumper.statbuf.st_size);
    //shc_log_info(shc_logname, "statbuf.st_size %ld\n", dumper.statbuf.st_size);
    uint64_t num_line = 0;
    uint64_t num_read = 0;
    char * header_ptr;
    char * seq_ptr;
    int seq_len;
    Block_timer read_prog_timer;
    start_timer(&read_prog_timer);

    {
        std::string message = "Assign "+ std::to_string(num_single_read_file)
                            + " single-ended read to components\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t progress_step = num_single_read_file/100;


        while(dumper.get_src_line(header_ptr, seq_ptr, seq_len))
        {
            assign_single_read_to_file_mmap(dumper ,header_ptr, seq_ptr, seq_len);
            num_read++;
            if((num_read-last_num) > progress_step)
            {
                int percentage = (100 * num_read)/num_single_read_file;
                progress.write(static_cast<double>(num_read) /
                               static_cast<double>(num_single_read_file));

                shc_log_info(shc_logname, "[ %u%] %u read out of %u is processed\n",
                        percentage, num_read, (num_single_read_file));
                log_stop_timer(&read_prog_timer);
                start_timer(&read_prog_timer);
                last_num = num_read;
            }
        }
    }
}


void Contig_graph_handler::assign_kmer_to_components()
{
    start_timer(&cgh_timer);
    shc_log_info(shc_logname, "Start assigning kmer to component\n");
    //creating dir
    typedef boost::filesystem::path path_t;
    typedef boost::filesystem::path dir_t;
    typedef boost::filesystem::path filename_t;
    char bases[33];
    bases[kmer_length] = '\0';

    add_or_overwrite_directory(lf.output_components_kmer_dir, lf.output_path);
    //create files
    kmer_dumper.setup_disk_dumper(lf.output_components_kmer_dir, lf.comp_kmer_prefix);

    uint64_t num_dumped_kmer = 0;
    uint64_t last_num_dumped = 0;
    uint64_t num_total_kmer = kh->kmer_counter.size();
    Block_timer kmer_prog_timer;
    start_timer(&kmer_prog_timer);
    Kmer_counter_map_iterator it;

    {
        std::string message = "Assign " + std::to_string(num_total_kmer)
                            + " kmer to components\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t progress_step = num_total_kmer/100;

        for(it=kh->kmer_counter.begin(); it!=kh->kmer_counter.end(); it++)
        {
            decode_kmer(bases, &(it->first), kmer_length);

            kmer_dumper.dump_kmer(component_array[it->second.contig], bases, it->second.count);
            //comp_num_kmers[component_array[it->second.contig]] ++;

            if(metis_setup.is_multiple_partition &&
            component_array_aux[it->second.contig] < curr_component_num)
            {
                kmer_dumper.dump_kmer(component_array_aux[it->second.contig], bases, it->second.count);
                //comp_num_kmers[component_array_aux[it->second.contig]] ++;
            }

            if((num_dumped_kmer -last_num_dumped) > progress_step)
            {
                int percentage = (100 * num_dumped_kmer)/num_total_kmer;
                progress.write(static_cast<double>(num_dumped_kmer) /
                               static_cast<double>(num_total_kmer));
                shc_log_info(shc_logname, "[ %u%] %u read out of %u is processed\n",
                        percentage, num_dumped_kmer, (num_total_kmer));
                log_stop_timer(&kmer_prog_timer);
                start_timer(&kmer_prog_timer);
                last_num_dumped = num_dumped_kmer;
            }
            num_dumped_kmer++;
        }
    }

    kmer_dumper.write_mem_to_disk();
    shc_log_info(shc_logname, "Finish assigning kmer to component\n");
}

void Contig_graph_handler::load_data_and_partition_reads()
{
    if(!exist_path(setting.local_files.output_comp_path))
    {
        shc_log_error("comp array data not available\n");
        exit(1);
    }

    std::ifstream fileReader(setting.local_files.output_comp_path);
    component_array.clear(); // index represents contig id, content stores comp id
    component_array_aux.clear(); // used for re-partition

    std::string num_comp_str, comp_str, comp_re_str, temp;
    int comp = 0;
    int comp_re = 0;

    std::getline(fileReader, num_comp_str, ' ');
    std::getline(fileReader, temp);
    curr_component_num = std::stoi(num_comp_str);

    shc_log_info(shc_logname, "num comp %d\n", curr_component_num);
    std::cout << "num comp " << curr_component_num << std::endl;


    if(!metis_setup.is_multiple_partition)
    {
        while(  std::getline(fileReader, comp_str, '\t') &&
                std::getline(fileReader, temp))
        {
            comp = std::stoi(comp_str);
            component_array.push_back(comp);
        }
    }
    else
    {
        while(  std::getline(fileReader, comp_str, '\t') &&
                std::getline(fileReader, comp_re_str , '\t') &&
                std::getline(fileReader, temp))
        {
            comp = std::stoi(comp_str);
            comp_re = std::stoi(comp_re_str);

            //shc_log_info(shc_logname, "ori %d, re %d\n", comp, comp_re);

            component_array.push_back(comp);
            component_array_aux.push_back(comp_re);
        }
    }
    /*
    shc_log_info(shc_logname, "comp array\n");
    for(std::vector<comp_num_t>::iterator it= component_array.begin();
                                   it!=component_array.end(); it++)
    {
        shc_log_info(shc_logname, "%d\n", *it);
    }

    shc_log_info(shc_logname, "aux comp array\n");
    for(std::vector<comp_num_t>::iterator it= component_array_aux.begin();
                                   it!=component_array_aux.end(); it++)
    {
        shc_log_info(shc_logname, "%d\n", *it);
    }
    */
}

/*
void Contig_graph_handler::log_contig_graph_edge(bool has_detail)
{
    shc_log_info(shc_logname, "num vertex: %d, num edge: %d \n",
                          boost::num_vertices(graph), boost::num_edges(graph));
    if(has_detail)
    {
        eip_t eip = boost::edges(graph);
        for(ei_t ei=eip.first; ei!=eip.second; ei++)
        {
            shc_log_info(shc_logname, "w: %d\n", graph[*ei].weight);
        }
    }
}

void Contig_graph_handler::log_comp_content()
{
    typedef std::map<comp_num_t, std::vector<contig_num_t> > comp_contig_map_t;
    typedef std::map<comp_num_t, std::vector<contig_num_t> >::iterator comp_contig_map_iterator;

    comp_contig_map_t comp_contig_map;
    for(int i=0; i<component_array.size(); i++)
    {
        comp_num_t comp_i = component_array[i];
        if(comp_i!=IMPOSSIBLE_COMP_NUM)
        {
            comp_contig_map[comp_i].push_back(i);
        }
    }

    for(comp_contig_map_iterator it = comp_contig_map.begin();
                                 it!=comp_contig_map.end(); ++it)
    {
        shc_log_info(shc_logname, "comp %d has %d contigs\n", it->first, it->second.size());
    }
}

void Contig_graph_handler::log_contigs(std::string file_path)
{
    std::ofstream file_writer(file_path);
    std::pair<vi_t, vi_t> vip = boost::vertices(graph);
    //walk through the graph list, and update if
    for(vi_t it=vip.first; it!=vip.second; it++)
    {
        bundled_contig_index & node = graph[*it];
        contig_num_t contig_i = node.contig_index;
        Contig_handler::Contig_base_info contig_info = ch->get_contig(contig_i);
        std::string contig(contig_info.base_start, contig_info.len);
        file_writer << (contig) << std::endl;
    }
    file_writer.close();
}
*/
