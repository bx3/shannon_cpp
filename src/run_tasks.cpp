#include "run_tasks.h"

Block_timer timer;

void set_stack_limit()
{
    const rlim_t kStackSize = 1024 * 1024 * 1024;   // min stack size = 16 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                std::cout << "setrlimit returned result = " << result << std::endl;
            }
        }
    }
}

void show_main_timer(struct Block_timer * main_timer)
{
    std::cout << "Since program starts, used ";
    stop_timer(main_timer);

    shc_log_info(shc_logname, "Since program starts, used\n");
    log_stop_timer(main_timer);
}

void load_default_setting(Shannon_C_setting & setting, std::string setting_path)
{
    boost::filesystem::path base_path = boost::filesystem::current_path();
    std::string base_path_str(base_path.c_str());

    parser_setting_file(setting_path, setting);
    memcpy(shc_logname, setting.local_files.log_filename_path.c_str(),
    setting.local_files.log_filename_path.size());

    if(!boost::filesystem::exists(setting_path))
    {
        shc_log_error("No default setting file: %s\n", setting_path.c_str());
        exit(1);
    }
}

void run_entire(int argc, char** argv)
{
    Shannon_C_setting setting;
    //load_default_setting(setting);



}

void run_custom(int argc, char** argv, Shannon_C_setting setting)
{
    int desiredTest = -1;
    while(desiredTest<0 || desiredTest>21)
    {
        std::cout << "Select desired application: \n";
        std::cout << "   0) specific test\n";
        std::cout << "   1) test_encoding_decoding\n";
        std::cout << "   2) test_Kmer with list\n";
        std::cout << "   3) test_Contig_graph\n";
        std::cout << "   4) test load and process contig_graph\n";
        std::cout << "   5) test_load_dump\n";
        std::cout << "   6) test pair read contig graph\n";
        std::cout << "   7) test load and process pair read contig graph\n";
        std::cout << "   8) test sequence graph\n";
        std::cout << "   9) test test_multi_seq_graph\n";
        std::cout << "   10) test sparse flow\n";
        std::cout << "   11) test multi-thread sparse flow\n";
        std::cout << "   12) test all\n";
        std::cout << "   13) test test_multi_seq_sf\n";
        std::cout << "   14) test test_Contig_graph with sorted kmer\n";
        std::cout << "   15) test_assigning_reads_and_kmers\n";
        std::cout << "   16) test_evaluation\n";
        std::cout << "   17) test_seq_run_it\n";
        std::cout << "   18) sort kmer\n";
        std::cout << "   19) test_filter_output\n";
        std::cout << "   20) test load_contig_patition_and_multi_seq_sf\n";
        std::cout << "   21) test sorted_kmer_contig_and_multi_seq_sf\n";
        std::cin >> desiredTest;
    }


    info_log_info(setting.local_files.timing_path.c_str(), "Command line Inputs\n");
    std::string cmd;
    for(int i=0; i<argc; i++)
    {
        std::string sub_cmd = argv[i];
        cmd += sub_cmd + " ";
    }
    info_log_info(setting.local_files.timing_path.c_str(), "%s\n\n", cmd.c_str());
    info_log_info(setting.local_files.timing_path.c_str(), "Desired test: %d\n\n", desiredTest);

    print_and_log_all_setting(setting);

    std::string log_suffix = "test" + desiredTest;
    setting.local_files.log_filename_path += log_suffix;

    shc_log_info(shc_logname, "Shannon C is starting\n");
    switch(desiredTest){
        case 0:
            test_specific(&setting);
            break;
        case 1:
            test_encoding_decoding();
            break;
        case 2:
            test_Kmer(&setting);
            break;
        case 3:
            test_Contig_graph(&setting);
            break;
        case 4:
            test_load_contig_graph(&setting);
            break;
        case 5:
            test_load_dump(&setting);
            break;
        case 6:
            test_pair_read_contig_graph(&setting);
            break;
        case 7:
            test_load_pair_contig_graph(&setting);
            break;
        case 8:
            test_seq_graph(&setting);
            break;
        case 9:
            test_multi_seq_graph(&setting);
            break;
        case 10:
            test_sparse_flow(&setting);
            break;
        case 11:
            test_multithread_sparse_flow(&setting);
            break;
        case 12:
            test_all(&setting);
            break;
        case 13:
            test_multi_seq_sf(&setting);
            break;
        case 14:
            test_sorted_kmer_contig(&setting);
            break;
        case 15:
            test_assigning_reads_and_kmers(&setting);
            break;
        case 16:
            test_evaluation(&setting);
            break;
        case 17:
            test_custom_seq_graph(&setting);
            break;
        case 18:
            sort_kmer(&setting);
            break;
        case 19:
            test_filter_output(&setting);
            break;
        case 20:
            test_load_contig_patition_and_multi_seq_sf(&setting);
            break;
        case 21:
            test_sorted_kmer_contig_and_multi_seq_sf(&setting);
            break;
        default:
            break;
    }

    return;
}

void test_sorted_kmer_contig_and_multi_seq_sf(Shannon_C_setting * setting)
{
    //get_num_seq(setting->local_files.reconstructed_seq_path);
    {
        test_sorted_kmer_contig(setting);
    }
    test_multi_seq_sf(setting);
    return;
}

void test_filter_output(Shannon_C_setting * setting)
{
    {
        Multi_graph_handler multi_graph(*setting);
        multi_graph.combine_all_reconst_seq_to_output();
        multi_graph.filter_output_seq_by_length();
        multi_graph.filter_reconstructed();
    }
    eval_reconstructed_seq(setting->local_files);
}

void sort_kmer(Shannon_C_setting * setting)
{
    start_timer(&timer);
    Duplicate_setting & dup = setting->dup_setting;
    Local_files * lf = &setting->local_files;

    Kmer_handler kmer_handler(*setting);

    Contig_handler contig_handler(setting->is_compress);
    kmer_handler.get_contig_handler(&contig_handler);
    kmer_handler.sort_kmer_descending_count_external();
}

void test_evaluation(Shannon_C_setting * setting)
{
    Multi_graph_handler multi_graph(*setting);
    eval_reconstructed_seq(setting->local_files);
}

void test_assigning_reads_and_kmers(Shannon_C_setting * setting)
{
    Duplicate_setting & dup = setting->dup_setting;
    Metis_setup & metis_setup = setting->metis_setup;
    Local_files * lf = &setting->local_files;

    shc_log_info(shc_logname, "Shannon C loads previous files\n");
    Kmer_handler kmer_handler(lf);
    Contig_handler contig_handler(setting->is_compress);
    kmer_handler.get_contig_handler(&contig_handler);

    kmer_handler.load_kmers_with_info_from_file(lf->output_kmer_path);

    Contig_graph_handler graph_handler(kmer_handler.get_kmer_length()-1,
                                       &kmer_handler,
                                       &contig_handler, *setting);
    graph_handler.load_data_and_partition_reads();
    Block_timer timer;
    start_timer(&timer);
    add_or_overwrite_directory(lf->output_components_read_dir, lf->output_path);
    if(setting->has_single)
        graph_handler.assign_reads_to_components_mmap(setting->local_files.input_read_path);
    if(setting->has_pair)
    {
        //std::cout << "process pair end " << std::endl;
        graph_handler.assign_paired_read_to_components_mmap(
                            setting->local_files.input_read_path_1,
                            setting->local_files.input_read_path_2);
    }

    //std::cout << "finish assign_reads_to_components with mmap" << std::endl;
    //stop_timer(&timer);
    add_or_overwrite_directory(lf->output_components_kmer_dir, lf->output_path);
    graph_handler.assign_kmer_to_components();
}



void test_sorted_kmer_contig(Shannon_C_setting * setting)
{
    start_timer(&timer);
    Duplicate_setting & dup = setting->dup_setting;
    Local_files * lf = &setting->local_files;

    Kmer_handler kmer_handler(*setting);

    Contig_handler contig_handler(setting->is_compress);
    kmer_handler.get_contig_handler(&contig_handler);

    kmer_handler.run_kmer_handler();

    kmer_handler.dump_kmers_with_info_to_file(lf->output_kmer_path);
    std::cout << "after dumping kmer" <<  std::endl;
    contig_handler.dump_all_contig(lf->output_contig_path);
    std::cout << "after dumping contig" <<  std::endl;

    kmer_handler.clear_kmer_map();
    kmer_handler.load_kmers_with_info_from_file(lf->output_kmer_path);

    Contig_graph_handler graph_handler(setting->kmer_length-1,
                        &kmer_handler, &contig_handler, *setting);

    graph_handler.run_contig_graph_handler();

    //kmer_handler.deallocate_kmer_map();
    show_main_timer(&timer);
}

void test_multi_seq_sf(Shannon_C_setting * setting)
{
    {
        Block_timer main_timer;
        start_timer(&main_timer);
        // starting multibridge and sparse flow
        Multi_graph_handler multi_graph(*setting);
        multi_graph.run_multi_seq_graph();

        std::cout << "multi seq graph finish, " << std::endl;
        show_main_timer(&main_timer);

        // sparse flow
        multi_graph.run_multi_sparse_flow(-1);
        std::cout << "multi sf finish, " << std::endl;
        show_main_timer(&main_timer);

        multi_graph.combine_all_reconst_seq_to_output();
        multi_graph.filter_output_seq_by_length();
        //multi_graph.filter_FP(setting->local_files.reconstructed_seq_path);
        multi_graph.filter_reconstructed();

        shc_log_info(shc_logname, "\t\t** entire process finish, \n");
        show_main_timer(&main_timer);
    }

    eval_reconstructed_seq(setting->local_files);
}

void test_all(Shannon_C_setting * setting)
{
#ifdef USE_MLOCK
    std::cout << "Using memory lock " << std::endl;
    if(mlockall(MCL_CURRENT||MCL_FUTURE)== -1)
    {
        perror("mlockall:1");
        exit(1);
    }
#endif


    Block_timer main_timer;
    start_timer(&main_timer);

    Block_timer part_timer;
    start_timer(&part_timer);
    run_jellyfish(*setting);
    std::cout << "Jellyfish task finishes, using time " << std::endl;
    stop_timer(&part_timer);
    shc_log_info(shc_logname, "finish run jellyfish\n");
    log_stop_timer(&part_timer);


    info_log_info(setting->local_files.timing_path.c_str(),
                        "jellyfish finish at %u sec\n", take_time(&part_timer));
    start_timer(&part_timer);

    Duplicate_setting & dup = setting->dup_setting;
    Local_files * lf = &setting->local_files;

    {
        // Getting graph be partitioned
        Kmer_handler kmer_handler(*setting);
        Contig_handler contig_handler(setting->is_compress);
        kmer_handler.get_contig_handler(&contig_handler);

        kmer_handler.run_kmer_handler();

        //std::cout << "start dump filtered kmer" << std::endl;
        Block_timer kmer_dump_timer;
        start_timer(&kmer_dump_timer);
        kmer_handler.dump_kmers_with_info_to_file(lf->output_kmer_path);
        shc_log_info(shc_logname, "finish dump filtered kmer\n");
        log_stop_timer(&kmer_dump_timer);
        //std::cout << "finish dump filtered kmer" << std::endl;
        //stop_timer(&kmer_dump_timer);

        //std::cout << "start dump contig" << std::endl;
        Block_timer contig_dump_timer;
        start_timer(&contig_dump_timer);
        contig_handler.dump_all_contig(lf->output_contig_path);
        shc_log_info(shc_logname, "finish dump contig\n");
        log_stop_timer(&contig_dump_timer);
        //std::cout << "finish dump contig" << std::endl;
        //stop_timer(&contig_dump_timer);
        //kmer_handler.clear_kmer_map();
    }

    info_log_info(setting->local_files.timing_path.c_str(),
                        "building contig finish at %u sec\n", take_time(&main_timer));
    start_timer(&part_timer);

    // reload everything just to free memory
    {
        show_main_timer(&main_timer);

    #ifdef USE_MLOCK
        if(munlockall()== -1)
        {
            perror("munlockall:1");
            exit(1);
        }
        std::cout << "Stop using memory lock " << std::endl;
    #endif

        //      create ad partition contig graph
        //kmer_handler.deallocate_kmer_map();

        Kmer_handler kmer_handler(lf);
        Contig_handler contig_handler(setting->is_compress);
        kmer_handler.get_contig_handler(&contig_handler);

        kmer_handler.load_kmers_with_info_from_file(lf->output_kmer_path);
        contig_handler.load_contig_file(lf->output_contig_path);

        Contig_graph_handler contig_graph_handler(setting->kmer_length-1,
                            &kmer_handler, &contig_handler, *setting);
        contig_graph_handler.run_contig_graph_handler();

        //kmer_handler.deallocate_kmer_map();

        show_main_timer(&main_timer);

    }

    info_log_info(setting->local_files.timing_path.c_str(),
                        "partition finish at %u sec\n", take_time(&main_timer));
    start_timer(&part_timer);

    // starting multibridge and sparse flow
    {
        Multi_graph_handler multi_graph(*setting);
        multi_graph.run_multi_seq_graph();

        info_log_info(setting->local_files.timing_path.c_str(),
                            "multi-graph multi-bridge finish at %u sec\n", take_time(&main_timer));
        show_main_timer(&main_timer);

        // sparse flow
        multi_graph.run_multi_sparse_flow(-1);
        info_log_info(setting->local_files.timing_path.c_str(),
                            "multi-graph sparse-flow finish at %u sec\n", take_time(&main_timer));
        start_timer(&part_timer);

        multi_graph.combine_all_reconst_seq_to_output();
        multi_graph.filter_output_seq_by_length();
        //multi_graph.filter_FP(setting->local_files.reconstructed_seq_path);
        multi_graph.filter_reconstructed();
        info_log_info(setting->local_files.timing_path.c_str(),
                            "find-rep finish at %u sec\n", take_time(&main_timer));
        start_timer(&part_timer);

        shc_log_info(shc_logname, "\t\t** entire process finish, \n");
        show_main_timer(&main_timer);
    }

    info_log_info(setting->local_files.timing_path.c_str(),
                        "shannon finish at %u sec\n", take_time(&main_timer));
    start_timer(&part_timer);
    eval_reconstructed_seq(setting->local_files);
}

void test_load_contig_patition_and_multi_seq_sf(Shannon_C_setting * setting)
{
    {
        test_load_contig_graph(setting);
    }
    test_multi_seq_sf(setting);
}

/*
uint64_t get_num_sam_line(std::string in_file)
{
    FILE *cmd;
    char result[1024];
    std::string grep_cmd("wc -l ");
    grep_cmd += in_file;
    cmd = popen(grep_cmd.c_str(), "r");
    if (cmd == NULL) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    while (fgets(result, sizeof(result), cmd)) {
        ;//printf("%s", result);
    }
    pclose(cmd);
    uint64_t num_seq = atoi(result);
    //std::cout << num_seq << std::endl;
    return num_seq;
}
*/
void test_specific(Shannon_C_setting * setting)
{
    //Multi_graph_handler multi_graph(*setting);
    /*
    std::string samfile_path = "/data1/bowen/Shannon_D_seq/output_Vip_Chat_entire/ana/out.sam"; //"/data1/bowen/Shannon_D_seq/output_L4_Arf5_entire/hista_out/out.sam";
    std::cout << "sam path " << samfile_path << std::endl;
    double gr_thresh = 0.9;

    GR_filter gr_filter;
    gr_filter.parse_sam_get_matches(samfile_path);

    std::string gr_reconstructed_seq_path =
        setting->local_files.reconstructed_seq_path + "_gr_before_filter";
    std::string cmd("mv " + setting->local_files.reconstructed_seq_path+ " "
                          + gr_reconstructed_seq_path );
    run_command(cmd, true);

    std::ifstream reader(gr_reconstructed_seq_path);
    std::ofstream writer(setting->local_files.reconstructed_seq_path);

    std::string name, read_base, header;

    while(  std::getline(reader, header) &&
            std::getline(reader, read_base))
    {
        name = header.substr(1);
        double num_match = gr_filter.get_num_match(name);
        double ratio = num_match / read_base.size();
        shc_log_info(shc_logname, "%s has len %d, num_match %d ratio %f\n", name.c_str(),
                                    read_base.size(), ratio, (uint32_t)num_match);
        if(ratio > gr_thresh)
        {
            writer << header << std::endl;
            writer << read_base << std::endl;
        }
        else if(ratio > 1)
        {
            shc_log_error("ratio greater than 1\n");
            exit(0);
        }
        else if(ratio < 0 )
        {
            shc_log_error("ratio smaller than 0\n");
            exit(0);
        }

    }
    reader.close();
    writer.close();

    eval_reconstructed_seq(setting->local_files);
    */






    //Multi_graph_handler multi_graph(*setting);
    //multi_graph.filter_FP(setting->local_files.reconstructed_seq_path);


    /*
    Duplicate_setting & dup = setting->dup_setting;
    Metis_setup & metis_setup = setting->metis_setup;
    Local_files * lf = &setting->local_files;

    // parameter setting
    int start_mem = get_proc_mem_value();

    Kmer_handler kmer_handler(*setting);

    Contig_handler contig_handler(setting->is_compress);
    kmer_handler.get_contig_handler(&contig_handler);

    struct JF_stats jf_stats;
    parse_jf_info(*setting, jf_stats);
    */
    //start_timer(&timer);
    //kmer_handler.sort_kmer_descending_count();
    /*
    std::cout << "size of kmer " << sizeof(Kmer_info) << std::endl;

    Duplicate_setting & dup = setting->dup_setting;
    Metis_setup & metis_setup = setting->metis_setup;
    Local_files * lf = &setting->local_files;

    shc_log_info(shc_logname, "Shannon C loads previous files\n");
    Kmer_handler kmer_handler(lf, setting->dup_setting.multiple_size);
    Contig_handler contig_handler(setting->is_compress);
    kmer_handler.get_contig_handler(&contig_handler);

    Contig_graph_handler graph_handler(setting->kmer_length-1,
                                       &kmer_handler,
                                       &contig_handler, *setting);
    graph_handler.load_data_and_partition_reads();
     */
    /*
    // test read
    struct Collect_reads coll_read_list(false, true,
                     101,
                     0,
                     0,
                     false,
                     *setting);
    Single_read_list & s_reads = coll_read_list.s_reads;

    std::map<std::string, int> inputs;

    Local_files & lf = setting->local_files;
    std::ifstream file_reader(lf.output_components_read_dir+"/comp0");
    std::string base, temp;
    read_num_t read_index;
    while(std::getline(file_reader, temp),
          std::getline(file_reader, base))
    {
        inputs.insert(std::make_pair(base,  1));
        s_reads.add_read(base, read_index);
    }

    for(read_num_t i=0; i<coll_read_list.get_num_reads(); i++)
    {
        Read_acc acc = s_reads.get_read(i);
        std::string read(acc.read_ptr, acc.len);
        if(inputs.find(read) == inputs.end())
        {
            std::cout << "see unknown read " << read << std::endl;
        }
    }

    std::cout << "num read " << s_reads.get_num_reads() << std::endl;
    std::cout << "total num read " << s_reads.get_total_num_read() << std::endl;
    */
}

void test_multithread_sparse_flow(Shannon_C_setting * setting)
{
    {
        int comp_i = -1;
        std::cout << "enter desired component to test, -1 for all component" << std::endl;
        std:: cin >> comp_i;

        Multi_graph_handler multi_graph(*setting);
        multi_graph.run_multi_sparse_flow(comp_i);

        multi_graph.combine_all_reconst_seq_to_output();
        multi_graph.filter_output_seq_by_length();
        //multi_graph.filter_FP(setting->local_files.reconstructed_seq_path);
        multi_graph.filter_reconstructed();
    }
    eval_reconstructed_seq(setting->local_files);
}

void test_sparse_flow(Shannon_C_setting * setting)
{
    int comp_i = 10;
    int graph_i = 0;
    std::cout << "enter desired component to test" << std::endl;
    std:: cin >> comp_i;
    std::cout << "enter desired graph to test" << std::endl;
    std:: cin >> graph_i;
    Sparse_flow_handler sparse_flow_handler(*setting);

    Local_files & lf = setting->local_files;

    std::string node_path(lf.output_seq_graph_path +
                                 lf.node_prefix + std::to_string(comp_i) +
                                 "/node" + std::to_string(graph_i));
    std::string edge_path(lf.output_seq_graph_path +
                                 lf.edge_prefix + std::to_string(comp_i) +
                                 "/edge" + std::to_string(graph_i));
    std::string path_path(lf.output_seq_graph_path +
                                 lf.path_prefix + std::to_string(comp_i) +
                                 "/path" + std::to_string(graph_i));
    std::string read_path(lf.output_seq_graph_path +
                          lf.read_prefix + std::to_string(comp_i) +
                          "/read" + std::to_string(graph_i));

    sparse_flow_handler.process_one_graph_component(comp_i, graph_i, node_path,
                        edge_path, path_path, read_path);
}



void test_multi_seq_graph(Shannon_C_setting * setting)
{
    if(boost::filesystem::exists(setting->local_files.reconstructed_sf_path))
        overwirte_a_file(setting->local_files.reconstructed_sf_path);
    Multi_graph_handler multi_graph(*setting);
    multi_graph.run_multi_seq_graph();
}

void test_custom_seq_graph(Shannon_C_setting * setting)
{
    uint64_t comp_i;
    std::cout << "enter desired component to test" << std::endl;
    std:: cin >> comp_i;

    Sequence_graph_handler seq_graph_handler(*setting);

}

void test_seq_graph(Shannon_C_setting * setting)
{
    uint64_t comp_i;
    std::cout << "enter desired component to test" << std::endl;
    std:: cin >> comp_i;
    Sequence_graph_handler seq_graph_handler(*setting);

    seq_graph_handler.run_it(comp_i, true);
    /*
    Sequence_graph_handler seq_graph_handler(*setting);
    seq_graph_handler.set_curr_comp(comp_i);
    seq_graph_handler.setup_input_file_path(comp_i);
    //seq_graph_handler.build_kmer_graph_from_reads();
    //std::cout << "finish build_kmer_graph_from_reads" << std::endl;

    seq_graph_handler.build_kmer_graph_from_edges();
    std::cout << "finish build_kmer_graph_from_edges" << std::endl;

    //shc_log_info(shc_logname, "before condense\n");
    //seq_graph_handler.log_term_array(false);


    shc_log_info(shc_logname,  "before condense num nodes is %d, num edges is %d\n",
            seq_graph_handler.get_num_nodes(), seq_graph_handler.get_num_edges());
    seq_graph_handler.condense_graph();
    shc_log_info(shc_logname,  "after condense num nodes is %d, num edges is %d\n",
            seq_graph_handler.get_num_nodes(), seq_graph_handler.get_num_edges());
    std::cout << "finish condense graph" << std::endl;

    shc_log_info(shc_logname, "log nodes\n");
    seq_graph_handler.log_all_node(true, true);
    shc_log_info(shc_logname,  "after condense num nodes is %d, num edges is %d\n",
            seq_graph_handler.get_num_nodes(), seq_graph_handler.get_num_edges());


    //shc_log_info(shc_logname, "after condense\n");
    //seq_graph_handler.log_all_node(true, false);

    //seq_graph_handler.remove_all_suspicious_nodes();
    //shc_log_info(shc_logname,  "%d nodes %d edges after destroying suspicious nodes\n",
    //        seq_graph_handler.get_num_nodes(), seq_graph_handler.get_num_edges());

    //std::cout << "finish remove suspicious" << std::endl;

    //seq_graph_handler.collapse_all();
    //shc_log_info(shc_logname,  "%d nodes %d edges after collapsing nodes\n",
    //        seq_graph_handler.get_num_nodes(), seq_graph_handler.get_num_edges());
    //std::cout << "finish collapse" << std::endl;

    //seq_graph_handler.log_edge_count();
    //shc_log_info(shc_logname, "QQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n");

    seq_graph_handler.update_xnode_set();
    std::cout << "after update xnode  " << std::endl;
    seq_graph_handler.assign_read_to_xnode();
    std::cout << "finish assign read to xnode  " << std::endl;

    seq_graph_handler.bridge_all_xnodes();
    std::cout << "finish bridge_all_xnodes, " <<(seq_graph_handler.get_num_nodes()) << " nodes remain" << std::endl;
    shc_log_info(shc_logname, "finish bridge_all_xnodes, %d node remain\n",
                                        seq_graph_handler.get_num_nodes());

    shc_log_info(shc_logname, "after multibridge\n");
    seq_graph_handler.log_all_node(true, false);
    shc_log_info(shc_logname, "finish bridge_all_xnodes, %d node remain\n",
                                        seq_graph_handler.get_num_nodes());

    //seq_graph_handler.find_approximate_copy_count();
    seq_graph_handler.break_self_loops();

    seq_graph_handler.break_all_cycles();


    //if(!seq_graph_handler.acyclic_check())
    //{
    //    shc_log_error("graph is not acyclic after break cycle\n");
    //    exit(1);
    //}

    //shc_log_info(shc_logname,"QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n");
    //std::cout << "number of nodes is " << seq_graph_handler.get_num_vertices() << std::endl;
    shc_log_info(shc_logname,  "after all num nodes is %d, num edges is %d\n",
            seq_graph_handler.get_num_nodes(), seq_graph_handler.get_num_edges());

    seq_graph_handler.find_known_path(setting->multi_graph_setup.max_hop_path);
    std::cout << "finish find all known path" << std::endl;

    std::cout << "finish all, " <<(seq_graph_handler.get_num_nodes()) << " nodes remain" << std::endl;
    shc_log_info(shc_logname, "finish all, %d node remain\n",
                                        seq_graph_handler.get_num_nodes());

    Local_files & lf = setting->local_files;

    std::string node_dir = lf.output_seq_graph_path + lf.node_prefix + std::to_string(comp_i);
    std::string edge_dir = lf.output_seq_graph_path + lf.edge_prefix + std::to_string(comp_i);
    std::string path_dir = lf.output_seq_graph_path + lf.path_prefix + std::to_string(comp_i);
    std::cout << "start output_components" << std::endl;
    if(boost::filesystem::exists(lf.reconstructed_seq_path))
        overwirte_a_file(lf.reconstructed_seq_path);
    seq_graph_handler.output_components(node_dir, edge_dir, path_dir);
    */
}

void test_load_pair_contig_graph(Shannon_C_setting * setting)
{
    shc_log_info(shc_logname, "Shannon C loads previous files\n");
    Kmer_handler kmer_handler(&setting->local_files);
    Contig_handler contig_handler(setting->is_compress);
    kmer_handler.get_contig_handler(&contig_handler);

    kmer_handler.load_kmers_with_info_from_file(setting->local_files.output_kmer_path);

    contig_handler.load_contig_file(setting->local_files.output_contig_path);

    Contig_graph_handler graph_handler(kmer_handler.get_kmer_length()-1,
                             &kmer_handler, &contig_handler,
                             *setting);

    graph_handler.group_components();


    graph_handler.dump_component_array(setting->local_files.output_comp_path);

    graph_handler.assign_paired_read_to_components(
                setting->local_files.input_read_path_1,
                setting->local_files.input_read_path_2);


    graph_handler.assign_kmer_to_components();
}

void test_pair_read_contig_graph(Shannon_C_setting * setting)
{
    Duplicate_setting & dup = setting->dup_setting;
    // parameter setting
    start_timer(&timer);
    Kmer_handler kmer_handler(*setting);

    Contig_handler contig_handler(setting->is_compress);

    //start processing
    kmer_handler.get_contig_handler(&contig_handler);
    std::cout << "Start reading kmer file"  << std::endl;
    if(kmer_handler.build_dict_from_kmer_file() != 0)
    {
        shc_log_error("Error reading kmers\n");
    }

    kmer_handler.sort_kmer_descending_count();

    kmer_handler.find_contig();

    kmer_handler.dump_kmers_with_info_to_file(setting->local_files.output_kmer_path);
    contig_handler.dump_all_contig(setting->local_files.output_contig_path);

    Contig_graph_handler graph_handler(kmer_handler.get_kmer_length()-1,
                             &kmer_handler, &contig_handler,
                              *setting);

    graph_handler.group_components();

    std::cout << "after group contig " <<  std::endl;


    if(setting->output_setup.read_to_component)
    graph_handler.assign_paired_read_to_components(
                setting->local_files.input_read_path_1,
                setting->local_files.input_read_path_2);
    if(setting->output_setup.kmer_to_component)
        graph_handler.assign_kmer_to_components();
    if(setting->output_setup.component_array)
        graph_handler.dump_component_array(setting->local_files.output_comp_path);

}

void test_load_dump(Shannon_C_setting * setting)
{
    Local_files * lf = &setting->local_files;
    print_and_log_local_file_system(lf);
    std::string kmer_back_file((lf->output_kmer_path) + "_re");
    std::string contig_back_file((lf->output_contig_path) + "_re");
    Kmer_handler kmer_handler(lf);
    kmer_handler.load_kmers_with_info_from_file(lf->output_kmer_path);
    kmer_handler.dump_kmers_with_info_to_file(kmer_back_file);

    Contig_handler contig_handler(setting->is_compress);
    contig_handler.load_contig_file(lf->output_contig_path);
    contig_handler.dump_all_contig(contig_back_file);
}

void test_load_contig_graph(Shannon_C_setting * setting)
{
    Duplicate_setting & dup = setting->dup_setting;
    Metis_setup & metis_setup = setting->metis_setup;
    Local_files * lf = &setting->local_files;

    shc_log_info(shc_logname, "Shannon C loads previous files\n");
    Kmer_handler kmer_handler(lf);
    Contig_handler contig_handler(setting->is_compress);
    kmer_handler.get_contig_handler(&contig_handler);

    kmer_handler.load_kmers_with_info_from_file(lf->output_kmer_path);

    contig_handler.load_contig_file(lf->output_contig_path);

    Contig_graph_handler graph_handler(kmer_handler.get_kmer_length()-1,
                                       &kmer_handler,
                                       &contig_handler, *setting);

    graph_handler.run_contig_graph_handler();
    /*
    graph_handler.group_components();


    graph_handler.dump_component_array(lf->output_comp_path);

    if(setting->output_setup.read_to_component)
        graph_handler.assign_reads_to_components(setting->local_files.input_read_path);

    if(setting->output_setup.kmer_to_component)
        graph_handler.assign_kmer_to_components();
    */
}

void test_Contig_graph(Shannon_C_setting * setting)
{
    Block_timer main_timer;
    start_timer(&main_timer);

    run_jellyfish(*setting);

    Duplicate_setting & dup = setting->dup_setting;
    Local_files * lf = &setting->local_files;

    {
        Kmer_handler kmer_handler(*setting);
        Contig_handler contig_handler(setting->is_compress);
        kmer_handler.get_contig_handler(&contig_handler);

        kmer_handler.run_kmer_handler();

        //std::cout << "start dump filtered kmer" << std::endl;
        Block_timer kmer_dump_timer;
        start_timer(&kmer_dump_timer);
        kmer_handler.dump_kmers_with_info_to_file(lf->output_kmer_path);
        shc_log_info(shc_logname, "finish dump filtered kmer\n");
        log_stop_timer(&kmer_dump_timer);
        //std::cout << "finish dump filtered kmer" << std::endl;
        //stop_timer(&kmer_dump_timer);

        //std::cout << "start dump contig" << std::endl;
        Block_timer contig_dump_timer;
        start_timer(&contig_dump_timer);
        contig_handler.dump_all_contig(lf->output_contig_path);
        shc_log_info(shc_logname, "finish dump contig\n");
        log_stop_timer(&contig_dump_timer);
        //std::cout << "finish dump contig" << std::endl;
        //stop_timer(&contig_dump_timer);
    }
    // Getting graph be partitioned

    {
        show_main_timer(&main_timer);

    #ifdef USE_MLOCK
        if(munlockall()== -1)
        {
            perror("munlockall:1");
            exit(1);
        }
        std::cout << "Stop using memory lock " << std::endl;
    #endif

        //      create ad partition contig graph
        //kmer_handler.deallocate_kmer_map();

        Kmer_handler kmer_handler(lf);
        Contig_handler contig_handler(setting->is_compress);
        kmer_handler.get_contig_handler(&contig_handler);

        kmer_handler.load_kmers_with_info_from_file(lf->output_kmer_path);
        contig_handler.load_contig_file(lf->output_contig_path);

        Contig_graph_handler contig_graph_handler(setting->kmer_length-1,
                            &kmer_handler, &contig_handler, *setting);
        contig_graph_handler.run_contig_graph_handler();

        //kmer_handler.deallocate_kmer_map();

        show_main_timer(&main_timer);

    }
}

void test_Kmer(Shannon_C_setting * setting)
{
    Duplicate_setting & dup = setting->dup_setting;
    Local_files * lf = &setting->local_files;

    // Getting graph be partitioned
    Kmer_handler kmer_handler(*setting);
    Contig_handler contig_handler(setting->is_compress);
    kmer_handler.get_contig_handler(&contig_handler);

    kmer_handler.run_kmer_handler();
}

/**
 *
 */
void test_encoding_decoding()
{
    /*
    uint8_t byte_list[4];
    char base_string[200] = "GGGCCAGGCACGGTGGCTCATGCCTATAATCCCAACACTTTGAGAGGCCAAGGCAGGTGGATCACT";
    char base_string_cp[200];
    char back[100];
    uint8_t length = 66;
    memcpy(base_string_cp, base_string, length);
    back[length] = '\0';
    encode_base_string(base_string, byte_list , length);
    //encode_in_place(base_string, (uint8_t*)base_string, length);
    decode_byte_list(byte_list , back, length);
    if(strcmp(back, base_string_cp) == 0 )
        printf("CORRECT, encode_base_string decode_base_string works \n");
    else
        printf("WRONG, encode_base_string decode_base_string\n");
        */
    char base_string[26] = "GGGCCAGGCACGGTGGCTCATGCCT";
    char base_string_cp[26];
    std::string reverse_str(base_string, 25);
    base_string_cp[25] = '\0';
    char back[26];
    back[25] = '\0';
    uint8_t length = 25;
    memcpy(base_string_cp, base_string, length);
    uint64_t byte;
    encode_reverse_kmer(base_string, &byte , length);
    std::reverse(reverse_str.begin(), reverse_str.end());
    uint64_t byte_s;
    encode_kmer(&reverse_str.at(0), &byte_s, length);
    if(byte_s==byte)
        std::cout << "correct" << std::endl;
    else
        std::cout << "wrong" << std::endl;

    decode_reverse_kmer(back, &byte, length);
    printf("decode %s\n", back);

    //uint64_t byte = 0;
    //encode_kmer(base_string, &byte, length);
    //print_bit(byte);
    //printf("start %s\n", base_string);
    //complement_num(&byte, length);
    //print_bit(byte);
/*
    decode_kmer( back, &byte, length);
    printf("back %s\n", back);

    if(combined == byte)
        printf("CORRECT, combining bytes equals to kmer direct encode\n");
    else
        printf("WRONG, combining bytes not equals to kmer direct encode\n");

    encode_kmer( base_string, &byte, length);
    decode_kmer( back, &byte, length);
    if(strcmp(back, base_string) == 0 )
        printf("CORRECT, encode/ decoder kmer directly works\n");
    else
        printf("WRONG, encode_base_string decode_base_string\n");
 * */
}
