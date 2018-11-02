#include "Multi_graph_handler.h"

int Multi_graph_handler::count_num_path_file(std::string path_dir)
{
    int num_path = 0;
    for(auto& entry : boost::make_iterator_range(
        boost::filesystem::directory_iterator(path_dir), {}))
    {
        std::string entry_str = entry.path().string();
        int underscore_index = 0;
        for(int i=entry_str.size()-1; i>=0; i--)
        {
            if(entry_str[i] == '_')
            {
                underscore_index = i;
                break;
            }
        }

        std::string path_prefix = "paths";
        std::string parse_str =
            entry_str.substr(underscore_index-path_prefix.size(), path_prefix.size());
        if(parse_str == path_prefix)
            num_path++;
    }
    return num_path;
}

int Multi_graph_handler::
count_num_component_sparse_flow()
{
    std::cout << "lf.output_seq_graph_path " << setting.local_files.output_seq_graph_path << std::endl;
    num_comp = count_num_path_file(setting.local_files.output_seq_graph_path);

    std::cout << "number of component for sparse flow is " << num_comp << std::endl;
    return num_comp;
}

int Multi_graph_handler::count_num_component_seq_graph()
{
    std::string kmer_path(setting.local_files.output_components_kmer_dir);
    int kmer_counter = count_num_files(kmer_path);
    //assert(read_counter==kmer_counter);
    num_comp = kmer_counter;
    std::cout << "number of component for multi-bridge is " << num_comp << std::endl;
    return num_comp;
}

void Multi_graph_handler::dump_all_single_nodes_to_reconstruct_output()
{
    std::string all_files_path(setting.local_files.single_node_dir + "/*");

    std::string cmd("find " + setting.local_files.single_node_dir +
                    " -name 'single_node_*' | sort -V | xargs cat"
                    " > " + setting.local_files.reconstructed_single_path);
    run_command(cmd, true);

    std::string rm_cmd("rm " + all_files_path);
    std::cout << "remove all temp single node seq" << std::endl;
    run_command(rm_cmd, true);
}

void Multi_graph_handler::filter_reconstructed()
{
    Block_timer lt;
    start_timer(&lt);
    std::string unfilter_output(setting.local_files.reconstructed_seq_path + "_unfiltered");
    if(exist_path(setting.local_files.reconstructed_seq_path))
    {
        std::string cmd("mv " + setting.local_files.reconstructed_seq_path + " "
                              + unfilter_output);
        run_command(cmd, true);
    }
    assert(exist_path(unfilter_output));


    find_representatives(unfilter_output, setting.local_files.reconstructed_seq_path);
    std::cout << "filter reconstructed fasta takes " <<std::endl;
    stop_timer(&lt);
    //deallocate_mem();
    shc_log_info(shc_logname, "finish filter_reconstructed\n");
    log_stop_timer(&lt);
}

void Multi_graph_handler::deallocate_mem()
{
    Rmer_contig_map().swap(rmer_contig_map);
}

void Multi_graph_handler::check_if_required_files_exist()
{
    if(exist_path(setting.local_files.output_components_read_dir))
        return;
    else
    {
        std::cerr << "[Error] partitioned components reads directory not exist at " << std::endl;
        std::cerr << setting.local_files.output_components_read_dir << std::endl;
        exit(0);
    }

    if(exist_path(setting.local_files.output_components_kmer_dir))
        return;
    else
    {
        std::cerr << "[Error] partitioned components kmer directory not exist at " << std::endl;
        std::cerr << setting.local_files.output_components_kmer_dir << std::endl;
        exit(0);
    }
}

void Multi_graph_handler::run_multi_seq_graph()
{
    struct Block_timer msq_timer;

    start_timer(&msq_timer);
    check_if_required_files_exist();

    int num_parallel = setting.num_parallel;
    int num_comp_for_seq = get_num_components_seq_graph();
    if(num_parallel > num_comp_for_seq)
        num_parallel = num_comp_for_seq;
    process_multi_seq_graph(num_parallel, num_comp_for_seq);
    dump_all_single_nodes_to_reconstruct_output();
    shc_log_info(shc_logname, "finish multibridge\n");
    log_stop_timer(&msq_timer);
#ifdef SHOW_PROGRESS
    std::cout << "finish multibridge, ";
    stop_timer(&msq_timer);
#endif
    //copy_file(setting.local_files.reconstructed_seq_path,
    //        setting.local_files.reconstruct_single_path);
}

// collect sf seq and single seq into final output
void Multi_graph_handler::collect_process_file(int num_parallel)
{
    overwirte_a_file(setting.local_files.reconstructed_sf_path);


    if(num_parallel > 0)
    {
        std::string cmd("find " + setting.local_files.output_seq_graph_result_path +
                    " -name '*.fasta'| sort -V | xargs cat");
        cmd += " > " + setting.local_files.reconstructed_sf_path;
        //std::cout << "cmd " << cmd << std::endl;
        run_command(cmd, true);

        //rm all temp process
        //for(int i=0; i<num_parallel; i++)
        //{
        //    std::string rm_cmd("rm ");
        //    std::string file_path =
        //                setting.local_files.output_seq_graph_result_path +
        //                "/process"+std::to_string(i);
        //    rm_cmd += file_path;
        //    run_command(rm_cmd, false);
        //}
        //std::cout << "remove all temp process sparse flow output file under comp_graph_output" << std::endl;
    }
}

void Multi_graph_handler::combine_contigs_seq_output()
{
    //std::string awk_cmd("awk '{if ( NR >= \"2\" ) print}' ");
    std::string cmd("cat " + setting.local_files.filtered_contig +
                    " >> " + setting.local_files.reconstructed_seq_path);
    run_command(cmd, true);
}

void Multi_graph_handler::combine_sf_single_seq_output()
{
    std::string seq_out_single(setting.local_files.reconstructed_single_path);
    if(boost::filesystem::exists(seq_out_single))
    {
        std::string cmd("cat  " + seq_out_single+ " " + setting.local_files.reconstructed_sf_path);
        cmd += " > ";
        cmd += setting.local_files.reconstructed_seq_path;
        run_command(cmd, true);
    }
}

void Multi_graph_handler::combine_all_reconst_seq_to_output()
{
    std::string single_seq = " ";
    std::string contig_seq = " ";
    //must have sparse flow output, even it is empty
    if(!boost::filesystem::exists(setting.local_files.reconstructed_sf_path))
    {
        std::cerr << "sparse flow output file [reconstructed_seq.fasta_sf] does not exist\n";
        exit(0);
    }

    // just leave
    if(!setting.take_single_node_seq && !setting.take_contig_seq)
    {
        std::string msg("Only take sparse flow output as final output");
        print_important_notice(msg);
        std::string cmd("cp " + setting.local_files.reconstructed_sf_path + " " +
                         setting.local_files.reconstructed_seq_path);
        run_command(cmd, true);
        return;
    }

    if(setting.take_single_node_seq)
    {
        single_seq = setting.local_files.reconstructed_single_path;
        if(!boost::filesystem::exists(single_seq))
        {
            std::cerr << "multibridge output single nodes file [reconstructed_seq.fasta_single] does not exist\n";
            exit(0);
        }
    }
    if(setting.take_contig_seq)
    {
        contig_seq = setting.local_files.filtered_contig;
        if(!boost::filesystem::exists(setting.local_files.filtered_contig))
        {
            std::cerr << "filtered contig output file [contig_filtered.fasta] does not exist\n";
            exit(0);
        }
    }

    std::string cmd("cat  " + single_seq + " " +
                     setting.local_files.reconstructed_sf_path + " " +
                     contig_seq);
    cmd += " > ";
    cmd += setting.local_files.reconstructed_seq_path;
    run_command(cmd, true);
}

// specific comp set to be -1 for all comp
void Multi_graph_handler::run_multi_sparse_flow(int specific_comp)
{
    struct Block_timer msf_timer;
    start_timer(&msf_timer);
    std::cout << "start multi-thread sparse flow" << std::endl;

    int num_parallel = setting.num_parallel;
    int num_comp_parse_flow = get_num_components_sparse_flow();
    process_sparse_flow_graph(num_parallel, num_comp_parse_flow, specific_comp);
    //exit(0);
    collect_process_file(num_parallel);

    shc_log_info(shc_logname, "finish sparse flow\n");
    log_stop_timer(&msf_timer);
#ifdef SHOW_PROGRESS
    std::cout << "finish sparse flow, ";
    stop_timer(&msf_timer);
#endif
}

void Multi_graph_handler::sort_work_list_by_size(std::deque<int> & work_list)
{
    std::vector<ID_size_pair> work_size_list;

    File_sorter file_sorter;
    std::sort(work_size_list.begin(), work_size_list.end(), file_sorter);
    for(int i=0; i<work_size_list.size(); i++)
    {
        work_list[i] = work_size_list[i].first;
    }
}

void Multi_graph_handler::sort_work_list_by_size(std::list<Comp_size_info> & work_list)
{
    std::ofstream writer(setting.local_files.component_mem_info);
    writer << "comp_id" << "\t\t" << "num_kmer" << "\t\t"
           << "read_size(byte)" << "\t\t"
           << "estimated_required_mem(byte)"
           << std::endl << std::endl;

    writer << "estimated_required_mem(byte) = (read size) + (num kmer)*"
           << std::to_string(ROUGH_EMPIRICAL_KMER_NUM_MEM_RATIO)  << std::endl
           << "since the last column is an estimate, the actual available "
           << "memory need to be " + std::to_string(mem_safe_ratio) << " times of the estimate "
           << "for algorithm to run without threat of crushing"
           << "either increase avail_mem arg or re-run partition algorithm "
           << "with smaller \"partition_size\" and \"partition_size\""
           << std::endl << std::endl;
    work_list.sort();
    for(std::list<Comp_size_info>::iterator it=work_list.begin();
                    it!=work_list.end(); it++)
    {
        int64_t estimated_mem = it->read_size + (it->kmer_num * ROUGH_EMPIRICAL_KMER_NUM_MEM_RATIO);
        writer << (it->comp_id) << "\t\t" << (it->kmer_num)
                  << "\t\t" << (it->read_size) << "\t\t" << estimated_mem << std::endl;
    }
    writer.close();
}

void Multi_graph_handler::
partition_work_to_process(int num_parallel, std::deque<int> & work_list,
                                std::vector<std::deque<int>> & process_queue)
{
    std::vector<int> place_order;
    for(int i=0; i<num_parallel; i++ )
        place_order.push_back(i);
    for(int i=num_parallel-1; i>=0; i-- )
        place_order.push_back(i);

    int j = 0;
    for(int i=0; i<work_list.size(); i++)
    {
        int work_i = work_list[i];
        int place_index = place_order[j];
        process_queue[place_index].push_back(work_i);
        if(j == place_order.size()-1)
            j = 0;
        else
            j++;
    }
}

bool Multi_graph_handler::get_comp_num_kmer(int i, int64_t & num)
{
    std::string file_path =
                     setting.local_files.output_components_kmer_dir
                     + setting.local_files.comp_kmer_prefix
                     + std::to_string(i);

    num = 0;
    if(!exist_path(file_path))
        return false;

    num = get_num_kmer(file_path);

    return true;
}

bool Multi_graph_handler::get_read_file_size(int i, int64_t & size)
{
    if(setting.has_single)
    {
        std::string file_path_prefix =
                         setting.local_files.output_components_read_dir
                         + setting.local_files.comp_read_prefix
                         + std::to_string(i);

        size = 0;
        int i=0;
        std::string file_path = file_path_prefix + "_" + std::to_string(i);
        if(!exist_path(file_path))
            return false;
        while(exist_path(file_path))
        {
            size += get_filesize(file_path);
            i++;
            file_path = file_path_prefix + "_" + std::to_string(i);
        }
        return true;
    }
    if(setting.has_pair)
    {
        std::string read_path_p1_prefix =
                     setting.local_files.output_components_read_dir
                     + setting.local_files.comp_read_prefix
                     + std::to_string(i)+"_p1";

        std::string read_path_p2_prefix =
                     setting.local_files.output_components_read_dir
                     + setting.local_files.comp_read_prefix
                     + std::to_string(i)+"_p2";
        if(!exist_path(read_path_p1_prefix + "_0") || !exist_path(read_path_p2_prefix+ "_0"))
            return false;

        size = 0;
        int i=0;
        std::string file_path_p1 = read_path_p1_prefix + "_" + std::to_string(i);
        while(exist_path(file_path_p1))
        {
            size += get_filesize(file_path_p1);
            i++;
            file_path_p1 = read_path_p1_prefix + "_" + std::to_string(i);
        }
        i = 0;
        std::string file_path_p2 = read_path_p2_prefix + "_" + std::to_string(i);
        while(exist_path(file_path_p2))
        {
            size += get_filesize(file_path_p2);
            i++;
            file_path_p2 = read_path_p2_prefix + "_" + std::to_string(i);
        }
        return true;
    }
}

void Multi_graph_handler::process_multi_seq_graph(int num_parallel, int num_components)
{
    Local_files & lf = setting.local_files;
    std::list<Comp_size_info> work_list;
    std::vector<pid_t> child_pids;
    int num, wk_fd, read_fd;

    int num_threads = 1;

    int start_class = 0;
    int end_class = num_components;
    for(int i=start_class; i<end_class; i++)
    {
        int64_t kmer_num, read_size;
        if(get_comp_num_kmer(i, kmer_num) &&
           get_read_file_size(i, read_size))
        {
            Comp_size_info p(i, kmer_num, read_size);
            work_list.push_back(p);
        }
        else
        {
            shc_log_error("%d file not exists\n", i);
            exit(1);
        }

    }
    sort_work_list_by_size(work_list);

    Comp_size_info largest_comp = work_list.front();
    if (!is_mem_enough_for_comp(largest_comp, ROUGH_EMPIRICAL_KMER_NUM_MEM_RATIO,
                                    avail_mem, mem_safe_ratio))
    {
        std::string message ("memory allocated " + std::to_string(setting.avail_mem)
                + " not enough for processing largest component, see file at \n"
                + setting.local_files.component_mem_info + "\n"
                + "to check component mem usage");
        print_important_notice(message);
        exit(0);
    }
    //exit(0);
    add_or_overwrite_directory(setting.local_files.output_seq_graph_path,
                               setting.local_files.output_path);
    std::cout <<  "setting.local_files.output_components_dump_read_dir " <<
                setting.local_files.output_components_dump_read_dir << std::endl;
    add_or_overwrite_directory(setting.local_files.output_components_dump_read_dir,
                               setting.local_files.output_path);

    //std::cout <<  "setting.local_files.output_components_thworn_read_dir " <<
    //       setting.local_files.output_components_thworn_read_dir << std::endl;
    //add_or_overwrite_directory(setting.local_files.output_components_thworn_read_dir,
    //                          setting.local_files.output_path);

    if(!setting.local_files.single_node_dir.empty())
    {
        //std::cout <<"single_node_dir " << setting.local_files.single_node_dir << std::endl;
        add_or_overwrite_directory(setting.local_files.single_node_dir,
                                       setting.local_files.output_path);
    }

    std::string wk_fifo = lf.fifo_dir + "/wk_fifo";
    std::vector<std::string> c_fifos;
    for(int i=0; i<num_parallel; i++)
    {
        std::string c_fifo = lf.fifo_dir + "/c_fifo" + std::to_string(i);
        c_fifos.push_back(c_fifo);
    }
    if (mkfifo(wk_fifo.c_str(), 0666) < 0)
        perror("mkfifo");
	for(int i=0; i<num_parallel; i++)
	{
	    if (mkfifo(c_fifos[i].c_str(), 0666) < 0)
   		    perror("mkfifo");
	}

    pid_t pid;
    std::cout << "parent pid = " << getpid() << std::endl;
    for (int i=0; i < num_parallel ; i++)
    {
        //std::cout << "create " << i << " children" << std::endl; // Child can keep going and fork once
        pid = fork();    // Fork
        //std::cout << "finish create " << i << " children" << std::endl; // Child can keep going and fork once
        if ( pid != 0) // if it is parent
        {
            child_pids.push_back(pid);
            usleep(5000);
            continue;
        }

        pid_t child_pid = getpid();
        std::cout << "Child pid " << child_pid << std::endl; // Child can keep going and fork once
        std::string child_mem_path = setting.local_files.mem_profiler.parallel_path_prefix
                                   + std::to_string(i);
        pid_t profiler_pid = fork_mem_profiler(child_pid, child_mem_path);
        // children process starts here
        int process_queue_index = i;

		if ((wk_fd = open(wk_fifo.c_str(), O_WRONLY)) < 0)
				perror("child - open");
		if ((read_fd = open(c_fifos[i].c_str(), O_RDONLY)) < 0)
				perror("child - open");

        pthread_mutex_t work_lock;
        if (pthread_mutex_init(&work_lock, NULL) != 0)
            { shc_log_error("unable to initialize mutex\n"); exit(1);}

        Seq_graph_works seq_graph_work;
        seq_graph_work.setting = setting;
        seq_graph_work.work_lock_ptr = &work_lock;
        seq_graph_work.process_i = process_queue_index;
        seq_graph_work.wk_fd = wk_fd;
        seq_graph_work.read_fd = read_fd;
        seq_graph_work.my_id  = std::to_string(process_queue_index);


        std::vector<pthread_t> threads(num_threads);

        for(int i=0; i<num_threads; i++)
        {
           //shc_log_info(shc_logname, "before create threads\n");
           int status = pthread_create(&threads.at(i), NULL,
               seq_graph_work.run_multi_seq_graph, (void*)(&(seq_graph_work)));
           usleep(5000);   // allow threads to be initialized

           if(status != 0)
           {
               std::cerr << "thread creation fails" << std::endl;
               exit(1);
           }
        }

        for(int i=0; i<threads.size(); i++)
        {
           pthread_join(threads[i], NULL);
        }
        close(wk_fd);
        close(read_fd);

        //pthread_mutex_destroy(&work_lock);

        //std::cout << "child process " << process_queue_index << " finish  " << std::endl;
        _exit(0);
    }

    usleep(500000);

    // parent
    if ((wk_fd = open(wk_fifo.c_str(), O_RDWR)) < 0)
        perror("parent - open");

    int c_fd[num_parallel];
    for(int i=0; i<num_parallel; i++)
    {
		if ((c_fd[i] = open(c_fifos[i].c_str(), O_WRONLY)) < 0)
        	perror("parent - open");
    }
    char write_buf[256];
    char read_buf[4096];
    int64_t num_process_finish = 0;
    int64_t num_comp_finish = 0;

    Process_comp_man proc_comp_man(num_parallel);

    std::string important_message("\nProgress bar below updates whenever one component finishes. "
                    "If there is one component with giant kmer-read files, it usually becomes "
                    "the bottleneck, and becomes the last component to finish. So if the process "
                    "bar stops at 95% or higher for long time, go to dir \n\n"
                    + setting.local_files.output_components_kmer_dir + "\n"
                    + setting.local_files.output_components_read_dir + "\n\n"
                    "check for file size (read files for each component are represented as a list 0,1,..) "
                    "to see if a few components are significantly greater\n\n"
                    + "Solution for slow processing time to a giant component:\n" +
                    + "1. change read sub-sampling parameter -g\n"
                    + "       for saving some multi-bridging attempts (preferred)\n"
                    + "2. reduce partition size: --metis_partition_size\n"
                    + "       for reducing component size, ( but might break long transcripts)\n");

    print_important_notice(important_message);

    {
        std::string message = "Multi-graph run " + std::to_string(work_list.size())
                            + " components\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t init_num_comp = work_list.size();
        std::vector<int> cli_list;

        do
        {
            if ((num = read(wk_fd, read_buf, sizeof(read_buf))) < 0)
    		    perror("parent - read");
            else
            {
                read_buf[num] = '\0';
    			std::string line_str(read_buf);
                int k = 0;

                std::vector<int> next_cli_list;

                for(int i=0; i<line_str.size(); i++)
                {
                    if (line_str[i] == ' ')
                    {
                        int len = i - k;
                        std::string pid_token = line_str.substr(k, len);
                        k = i + 1;

                        //std::cout << "pid_token " <<  pid_token << std::endl;

                        memcpy(write_buf, line_str.c_str() + k, len);
                        write_buf[len] = '\0';
                        int previous_comp = DEFAULT_COMP_NUM;
                        int delimitor_index = -1;
                        for(int j=0; j<pid_token.size(); j++)
                        {
                            if (pid_token[j] == ',')
                            {
                                delimitor_index = j;
                                break;
                            }
                        }
                        assert(delimitor_index!=-1);
                        std::string pid_str = pid_token.substr(0, delimitor_index);
                        std::string prev_mem_str = pid_token.substr(delimitor_index+1);

                        int cli_idx = boost::lexical_cast<int>(pid_str);
                        cli_list.push_back(cli_idx);

                        int64_t prev_mem = boost::lexical_cast<int64_t>(prev_mem_str);


                        if(prev_mem != DEFAULT_COMP_MEM)
                        {
                            num_comp_finish += 1;
                            Comp_size_info completed_comp = proc_comp_man.process_comp[cli_idx];
                            //std::cout << "comp " <<  completed_comp.comp_id << " prev_mem " << prev_mem << std::endl;
                            int64_t kmer_part = prev_mem;// - completed_comp.read_size ;
                            assert(kmer_part > 0 );
                            double kmer_num_mem_ratio = static_cast<double>(kmer_part) /
                                                 static_cast<double>(completed_comp.kmer_num);
                            proc_comp_man.kmer_num_mem_ratio_list.push_back(kmer_num_mem_ratio);
                            avail_mem += completed_comp.mem_required;
                            //std::cout << "avail_mem release " << avail_mem << std::endl;
                            //std::cout << "kmer_num_mem_ratio " << kmer_num_mem_ratio << std::endl;
                        }
    				}
                }


                double percentage = static_cast<double>(num_comp_finish) /
                                    static_cast<double>(init_num_comp);
                progress.write(percentage);

                //std::cout << "freshed cli "<< std::endl;
                //for(int i=0; i<cli_list.size(); i++)
                //{
                //    int fresh_cli = cli_list[i];
                //    std::cout << "process " << child_pids[fresh_cli] << " - " << fresh_cli << std::endl;
                //}


                for(int i=0; i<cli_list.size(); i++)
                {
                    double mean_ratio = proc_comp_man.get_mean_ratio();
                    //std::cout << "mean_ratio " << mean_ratio << std::endl;
                    Comp_size_info comp_info = get_next_comp(child_pids,
                                work_list, num_parallel, mean_ratio);

                    // when mem is not enough, not send anything to child and children
                    // chile block, save their request in the futurn and in
                    // higher priority
                    if(comp_info.comp_id == MEM_NOT_ENOUGH_ID)
                    {
                        next_cli_list.assign(cli_list.begin()+i, cli_list.end());
                        //next_cli_list.push_back(cli_list[i]);
                        std::cout << std::endl;
                        std::string message("\ndetect insufficient memory for " +
                                std::to_string(next_cli_list.size()) +
                                " un-processed component, child process(es) blocking, wait for memory release from any finishing child\n\n");
                        std::string blocked_pids(
                                        std::to_string(next_cli_list.size()) +
                                        " blocked process ");
                        for(int ni=0; ni<next_cli_list.size(); ni++)
                        {
                            int block_cli = next_cli_list[ni];
                            blocked_pids += std::to_string(child_pids[block_cli]) + " ";
                            //std::cout << "blocked process " <<  << " - " << block_cli << std::endl;
                        }
                        std::string print_message = message + blocked_pids + "\n\n";
                        print_important_notice(print_message);
                        break;

                    }

                    int cli_idx = cli_list[i];
                    int comp_to_process = comp_info.comp_id;
                    //std::cout << "comp_to_process " <<  comp_to_process << std::endl;

                    proc_comp_man.process_comp[cli_idx] = comp_info;

                    //process_comp[cli_idx].emplace_back(cli_idx, comp_info);

                    if (comp_to_process == -1)
                    {
                        num_process_finish++;
                    }
                    std::string comp_to_process_str = std::to_string(comp_to_process);

                    memcpy(write_buf, comp_to_process_str.c_str(), comp_to_process_str.size());
                    write_buf[comp_to_process_str.size()] = '\0';


                    if ((num = write(c_fd[cli_idx], write_buf, strlen(write_buf))) < 0)
                         perror("parent - write");
                }
                cli_list = next_cli_list;
            }


            if(num_process_finish == num_parallel*num_threads)
            {
                //printf("everything finishes\n");
                break;
            }

        } while (num > 0);
    }

    // collect zombie and last child process
    std::cout << "parent process is waiting for children " << std::endl;
    for(int i=0; i<num_parallel; i++)
    {
        int status;
        wait(&status);
    }
    std::cout << "parent process finish wait " << std::endl;

    close(wk_fd);

	for(int i=0; i<num_parallel; i++)
	{
		close(c_fd[i]);
   		unlink(c_fifos[i].c_str());
	}
	unlink(wk_fifo.c_str());
    return;
}

//return true when allow
bool Multi_graph_handler::
is_mem_enough_for_comp(Comp_size_info & curr_comp, double mean_ratio_,
                        int64_t avail_memory, double memory_safe_ratio)
{
    curr_comp.mem_required = mean_ratio_ * curr_comp.kmer_num + curr_comp.read_size;

    return curr_comp.mem_required*memory_safe_ratio < avail_memory;
}

Multi_graph_handler::Comp_size_info
Multi_graph_handler::get_next_comp(std::vector<pid_t> & child_pids,
        std::list<Comp_size_info> & work_list, int num_parallel, double mean_ratio)
{
    if(work_list.empty())
    {
        Comp_size_info a(-1,0,0);
        return a;
    }

    for(std::list<Comp_size_info>::iterator it=work_list.begin();
                        it!=work_list.end(); it++)
    {
        Comp_size_info curr_comp = *it;
        if(is_mem_enough_for_comp(curr_comp, mean_ratio, avail_mem, mem_safe_ratio))
        {
            //std::cout << "avail_mem       " << avail_mem << std::endl;
            avail_mem = avail_mem - curr_comp.mem_required;
            //std::cout << "mem_required    " << curr_comp.mem_required << std::endl;
            //std::cout << "avail_mem left  " << avail_mem << std::endl;
            work_list.erase(it);
            return curr_comp;
        }
    }

    // memory not enough for any comp
    Comp_size_info b(MEM_NOT_ENOUGH_ID,0,0);
    return b;
}

uint64_t Multi_graph_handler::get_parallel_total_mem(std::vector<pid_t> & child_pids, int num_parallel)
{
    uint64_t mem_sum = 0;
    for(int i=0; i<num_parallel; i++)
    {
        mem_sum += get_mem(child_pids[i]);
    }
    return mem_sum;
}

void Multi_graph_handler::process_sparse_flow_graph(int & num_parallel, int num_components, int specific_comp)
{
    int start_class = 0;
    int end_class = num_components;
    Local_files & lf = setting.local_files;
    if(specific_comp!=-1)
    {
        std::string comp_node_dir = lf.output_seq_graph_path + lf.node_prefix + std::to_string(specific_comp);
        std::string comp_edge_dir = lf.output_seq_graph_path + lf.edge_prefix + std::to_string(specific_comp);
        std::string comp_path_dir = lf.output_seq_graph_path + lf.path_prefix + std::to_string(specific_comp);
        if (exist_path(comp_node_dir) && exist_path(comp_edge_dir) &&
            exist_path(comp_path_dir))
        {
            start_class = specific_comp;
            end_class = specific_comp+1;
        }
        else
        {
            shc_log_error("please specify a valid component for multi-sparse flow\n");
            exit(1);
        }
    }

    std::cout << "start_class " << start_class <<std::endl;
    std::cout << "end_class " << end_class <<std::endl;
    std::cout << "num paralle " << num_parallel << std::endl;

    std::deque<Comp_graph> work_queue;

    for(int i=start_class; i<end_class; i++)
    {
        std::string path(lf.output_seq_graph_path +
                                     lf.edge_prefix + std::to_string(i));
        int num_sub_graph = count_num_files(path);
        for(int j=0; j<num_sub_graph; j++)
        {
            work_queue.push_back(Comp_graph(i, j));
        }
    }

    std::cout << "work_queue " << work_queue.size() << std::endl;

    if(num_parallel > work_queue.size())
        num_parallel = work_queue.size();

    std::vector<std::deque<Comp_graph>> process_queue(num_parallel, std::deque<Comp_graph>());

    partition_work_to_process_randomize(num_parallel, work_queue, process_queue);
    for(int i=0; i<process_queue.size(); i++)
    {
        //std::cout << "process " << i << ": ";
        std::deque<Comp_graph> & my_works = process_queue[i];
        std::cout << "process " << i << " has " << my_works.size()
                  << " components to process" << std::endl;;
        info_log_info(shc_logname, "process %d: ", i);
        for(std::deque<Comp_graph>::iterator it=my_works.begin();
                            it!=my_works.end(); it++)
        {
        //    std::cout << " ("<< (it->comp_i) << " " << it->graph_i << ") ";
            info_log_info(shc_logname, " (%d %d) ", (it->comp_i), it->graph_i);
        }
        //std::cout << std::endl;
        info_log_info(shc_logname, "\n");
    }

    //exit(0);
    add_or_overwrite_directory(setting.local_files.output_seq_graph_result_path,
                               setting.local_files.output_path);



    pid_t pid;
    int i;
    std::cout << "num parallel " << num_parallel << std::endl;
    std::cout << "parent pid = " << getpid() << std::endl;
    for (i=0; i < num_parallel ; i++)
    {
        //std::cout << "create " << i << " children" << std::endl; // Child can keep going and fork once
        pid = fork();    // Fork
        //std::cout << "finish create " << i << " children" << std::endl; // Child can keep going and fork once
        if ( pid != 0) // if it is parent
            continue;

        std::cout << "Child " << getpid() << std::endl; // Child can keep going and fork once

        // children process starts here
        int process_queue_index = i;
        std::deque<Comp_graph> my_works = process_queue[process_queue_index];
        pthread_mutex_t work_lock;
        if (pthread_mutex_init(&work_lock, NULL) != 0)
            { shc_log_error("unable to initialize mutex\n"); exit(1);}
        pthread_mutex_t write_lock;
        if (pthread_mutex_init(&write_lock, NULL) != 0)
            { shc_log_error("unable to initialize mutex\n"); exit(1);}


        Sparse_flow_works sparse_flow_work;
        sparse_flow_work.work_list = &my_works;
        sparse_flow_work.setting = setting;
        sparse_flow_work.work_lock_ptr = &work_lock;
        sparse_flow_work.write_lock_ptr = &write_lock;
        sparse_flow_work.process_i = process_queue_index;
        sparse_flow_work.output_path = setting.local_files.output_seq_graph_result_path +
                    "/process"+std::to_string(process_queue_index);
        sparse_flow_work.init_total_work = my_works.size();

        int num_threads = 1;
        std::vector<pthread_t> threads(num_threads);

        for(int i=0; i<num_threads; i++)
        {
           //shc_log_info(shc_logname, "before create threads\n");
           int status = pthread_create(&threads.at(i), NULL,
               sparse_flow_work.run_multi_sparse_flow, (void*)(&(sparse_flow_work)));
           usleep(5000);   // allow threads to be initialized

           if(status != 0)
           {
               std::cerr << "thread creation fails" << std::endl;
               exit(1);
           }
        }

        for(int i=0; i<threads.size(); i++)
        {
           pthread_join(threads[i], NULL);
        }
        pthread_mutex_destroy(&work_lock);
        pthread_mutex_destroy(&write_lock);
        std::cout << "child process " << process_queue_index << " finish  " << std::endl;
        _exit(0);
    }

    std::cout << "parent process is waiting for children " << std::endl;
    for(int i=0; i<num_parallel; i++)
    {
        int status;
        wait(&status);
    }
    std::cout << "parent process gets children " << std::endl;

    return;

}

void Multi_graph_handler::
partition_work_to_process_randomize(int num_process, std::deque<Comp_graph> & work_queue,
                            std::vector<std::deque<Comp_graph>> & process_queue)
{
    //std::cout << "before " << work_queue.size() << std::endl;
    auto engine = std::default_random_engine{};
    std::shuffle(std::begin(work_queue), std::end(work_queue), engine);
    //std::cout << "after  " << work_queue.size() << std::endl;

    std::vector<int> place_order;
    for(int i=0; i<num_process; i++ )
        place_order.push_back(i);
    for(int i=num_process-1; i>=0; i-- )
        place_order.push_back(i);

    int j = 0;
    for(int i=0; i<work_queue.size(); i++)
    {
        Comp_graph work_i = work_queue[i];
        int place_index = place_order[j];
        process_queue[place_index].push_back(work_i);
        if(j == place_order.size()-1)
            j = 0;
        else
            j++;
    }
}


bool Multi_graph_handler::
output_fasta_file_validator(std::string & output_path)
{
    std::ifstream reader(output_path);
    std::cout << "output path: " << output_path << std::endl;
    std::string line;

    int i = 0;

    if(system(NULL))
        ;
    else
        exit(EXIT_FAILURE);

    std::string cmd = "wc -l " + setting.local_files.reconstructed_seq_path;
    system(cmd.c_str());

    while(std::getline(reader, line))
    {
        if(line.empty())
            break;

        if(i%2 == 0 )
        {
            if(line.at(0) != '>')
            {
                shc_log_error("file has header without starting > \n");
                exit(1);
            }
        }
        else
        {
            if(line.at(0) == '>')
            {
                shc_log_error("file has seq with starting > \n");
                exit(1);
            }
        }

        i++;
    }

    if(i%2 != 0)
    {
        shc_log_error("file has uneven number of lines\n");
        exit(1);
    }
    else
    {
        std::cout << "it is a valid fasta file with " << i<< " lines" << std::endl;
        reader.close();
        return true;
    }

}


void Multi_graph_handler::
find_representatives(std::string in_file, std::string output_file)
{

    uint64_t num_seq;

    //std::cout << "getpid " << getpid() << std::endl;
    Block_timer timer;
    start_timer(&timer);

    if(exist_path(in_file))
        num_seq = get_num_seq(in_file);
    else
        num_seq = get_num_seq(setting.local_files.reconstructed_seq_path);

    int ave_read_length = 100;
    if(setting.has_single && !setting.has_pair)
        ave_read_length = setting.single_read_length;
    else if(!setting.has_single && setting.has_pair)
        ave_read_length =
              (setting.pair_1_read_length + setting.pair_2_read_length)/2;
    else if(setting.has_single && setting.has_pair)
        ave_read_length = (setting.pair_1_read_length +
                setting.pair_2_read_length + setting.single_read_length)/3;

    //std::vector<std::string> headers(num_seq, std::string(30, '\0'));
    //std::vector<std::string> seqs(num_seq, std::string(ave_read_length, '\0'));
    std::vector<std::string> headers;
    std::vector<std::string> seqs;

    std::string header, seq, line, rmer;
    rmer.resize(rmer_length);
    uint64_t contig_num = 0;
    rmer_contig_map.reserve(num_seq*ave_read_length);
    uint64_t index = 0;



    uint64_t byte;
    {
        std::ifstream file_reader(in_file);
        Seqs_header_map seqs_header_map;
        //seqs_header_map.reserve(num_seq);

        while( std::getline(file_reader, line))
        {
            if( line[0] == '>')
            {
                header = line.substr(1);
                //header = "hellio" + std::to_string(index++);
                headers.push_back(header);
                if(contig_num %SHOW_STEP == 0)
                {
                    std::cout << "processed " << contig_num << " sequences of out "
                              << num_seq <<  ", ";
                    stop_timer(&timer);
                    start_timer(&timer);
                }
            }
            else
            {
                if(seqs_header_map.find(line) != seqs_header_map.end())
                {
                    //Seqs_header_map_iterator it = ;
                    //std::cout << "find repeat seq " << std::endl;
                    //std::cout << seqs_header_map[line] << std::endl;
                    //std::cout << headers.back() << std::endl;

                    //std::cout << line << std::endl;
                    headers.pop_back();
                    continue;
                }

                if(line.size() > setting.output_seq_min_len)
                {
                    //std::cout << "hello0" << std::endl;
                    //std::cout << "hello1" << line.size() << std::endl;
                    seqs.push_back(line);
                    //std::cout << "hello2" << line.size() << std::endl;
                    for(uint64_t i=0; i<line.size()-rmer_length+1; i++)
                    {
                        encode_kmer(&line.at(i), &byte, rmer_length);
                        //std::cout << "hello3" << std::endl;
                        Seq_info seq_info(contig_num, i);
                        Rmer_contig_map_iterator it = rmer_contig_map.find(byte);
                        if(it != rmer_contig_map.end())
                        {
                            (it->second).push_back(seq_info);
                        }
                        else
                        {
                            std::vector<Seq_info> list_seq(1, seq_info);
                            list_seq.shrink_to_fit();//reserve(VEC_INIT_SIZE);
                            rmer_contig_map.insert(std::make_pair(byte, list_seq));
                        }
                    }
                    contig_num++;
                    seqs_header_map.insert(std::make_pair(line, headers.back()));
                }
                else
                {
                    headers.pop_back();
                }
            }
        }
        file_reader.close();
    }

    //std::cout << "after build rmer_contig_map " << std::endl;

    contig_num = 0;
    std::ofstream file_writer(output_file);
    Block_timer write_timer;
    start_timer(&write_timer);
    for(uint64_t i=0; i<seqs.size(); i++)
    {
        bool duplicate_suspect = duplicate_check_ends(seqs, i, false);
        if (setting.is_double_stranded && !duplicate_suspect)
        {
            duplicate_suspect = duplicate_check_ends(seqs, i, true);
        }
        if(contig_num++ %SHOW_STEP == 0)
        {
            int percentage = (100 * i)/seqs.size();
            std::cout << "[" << percentage << "%] "
                      << "written " << contig_num << " sequences from checked "
                      << i << "sequences out of " << seqs.size() << std::endl;
            stop_timer(&write_timer);
            start_timer(&write_timer);
        }
        if(!duplicate_suspect)
        {
            //std::cout << "header_index " << (header_index) << std::endl;
            file_writer << ">" + headers[i] << std::endl;
            file_writer << (seqs[i]) << std::endl;
        }
    }
    file_writer.close();
}

void Multi_graph_handler::filter_output_seq_by_length()
{
    std::string mv_cmd("mv ");
    mv_cmd = mv_cmd + setting.local_files.reconstructed_seq_path + " " +
                setting.local_files.unfiltered_length_reconst_seq;
    run_command(mv_cmd, true);
    std::ofstream filtered_writer(setting.local_files.reconstructed_seq_path);
    std::ifstream unfiltered_reader(setting.local_files.unfiltered_length_reconst_seq);
    std::string header, seq;
    while(std::getline(unfiltered_reader, header, '\n'),
          std::getline(unfiltered_reader, seq, '\n'))
    {
        if (seq.size() > setting.output_seq_min_len)
        {
            filtered_writer << header << std::endl;
            filtered_writer << seq << std::endl;
        }
    }
    filtered_writer.close();
    unfiltered_reader.close();
}

// return true if it is duplicate
bool Multi_graph_handler::
duplicate_check_ends(std::vector<std::string> & seqs, uint64_t header, bool rc)
{
    typedef std::pair<int64_t, int64_t> Index_pair;
    std::string seq = seqs[header];
    if (rc)
    {
        reverse_complement(seq);
    }
    uint64_t first_byte, last_byte;
    encode_kmer(&seq.at(0), &first_byte, rmer_length);
    encode_kmer(&seq.at(seq.size()-rmer_length), &last_byte, rmer_length);
    //std::string first_rmer(seq.substr(0, rmer_length));
    //std::string last_rmer(seq.substr(seq.size()-rmer_length, rmer_length));

    Rmer_contig_map_iterator first_rc_it = rmer_contig_map.find(first_byte);
    if(first_rc_it == rmer_contig_map.end())
        return false;

    Rmer_contig_map_iterator last_rc_it = rmer_contig_map.find(last_byte);
    if(last_rc_it == rmer_contig_map.end())
        return false;

    std::vector<Seq_info> & first_contig_vec = first_rc_it->second;
    std::vector<Seq_info> & last_contig_vec = last_rc_it->second;

    typedef tsl::hopscotch_map<uint64_t, Index_pair,
                   hash_u64, equ64> Contig_dict_map;
    typedef tsl::hopscotch_map<uint64_t, Index_pair,
                hash_u64, equ64>::iterator Contig_dict_map_iterator;
    Contig_dict_map contig_dict;

    for(std::vector<Seq_info>::iterator it=first_contig_vec.begin();
                            it!=first_contig_vec.end(); it++)
    {
        Seq_info & info = *it;
        if(info.header == header)
            continue;
        Contig_dict_map_iterator d_it = contig_dict.find(info.header);
        if(d_it != contig_dict.end())
        {
            (d_it.value()).first = info.index;
        }
        else
        {
            Index_pair index_pair(info.index, -1);
            contig_dict.insert(std::make_pair(info.header, index_pair));
        }
    }

    for(std::vector<Seq_info>::iterator it=last_contig_vec.begin();
                            it!=last_contig_vec.end(); it++)
    {
        Seq_info & info = *it;
        if(info.header == header)
            continue;
        Contig_dict_map_iterator d_it = contig_dict.find(info.header);
        if(d_it != contig_dict.end())
        {
            (d_it.value()).second = info.index;
        }
        else
        {
            Index_pair index_pair(-1, info.index);
            contig_dict.insert(std::make_pair(info.header, index_pair));
        }
    }

    for(Contig_dict_map_iterator it=contig_dict.begin();
                                it!=contig_dict.end(); it++)
    {
        Index_pair & index_pair = it.value();
        if(index_pair.first>=0 && index_pair.second>=0)
        {
            int64_t diff = index_pair.second - index_pair.first;
            int64_t seq_size = static_cast<int64_t>(seq.size());
            if(std::abs(diff+rmer_length-seq_size) < 3)
            {
                std::string & seq_other = seqs[it->first];
                if(seq.size()<seq_other.size() ||
                   (seq.size()==seq_other.size() && header > (it->first)) )
                {
                    return true;
                }
            }
        }
    }
    return false;
}


//take rec_fasta, read_1, read_2 files of type.
//Output in out_fasta. Temp dir = out_dir.
//Output fasta is in out_dir/reconstructed.fasta
//Flags is '-f --ff' for fasta and forward-forward
//Flags is '-q --fr' for fastq and forward-reverse
void Multi_graph_handler::filter_FP(std::string reconstructed_seq_path)
{
    /*
    if(setting.has_pair && setting.filter_paired)
    {
        shc_log_info(shc_logname, "start filter_FP\n");
        Local_files & lf = setting.local_files;
        std::string fp_dir = setting.local_files.output_path + "/fp_files";
        boost::filesystem::path fp_dir_boost(fp_dir);
        if( !boost::filesystem::exists(fp_dir_boost) )
            add_directory(fp_dir_boost);

        std::string output_hisat = fp_dir + "/rec.hisat";
        std::string output_sam = fp_dir + "/rec.sam";
        std::string output_bam = fp_dir + "/rec.bam";
        std::string output_sam_sort = fp_dir + "/rec_sort";
        std::string output_sam_depth = fp_dir + "/rec.depth";

        std::string & read_1 = lf.input_read_path_1;
        std::string & read_2 = lf.input_read_path_2;

        std::string hisat_build_cmd =
                "hisat-build " + reconstructed_seq_path + " " + output_hisat;
        std::cout << "Build hisat file" << std::endl;
        run_command(hisat_build_cmd, true);

        std::string hisat_align_cmd =
                "hisat --no-spliced-alignment --no-discordant -f --ff  -x "
                + output_hisat + " -1 " + read_1 + " -2 " + read_2 + " -S "
                + output_sam;
        std::cout << "Align" << std::endl;
        run_command(hisat_align_cmd, true);

        std::cout << "Process SAM / BAM file to get depth information" << std::endl;
        std::string samtool_bam_cmd = "samtools view -bS -f 0x2 " + output_sam +
                                      " > " + output_bam;
        run_command(samtool_bam_cmd, true);

        std::string sam_sort_cmd = "samtools sort " + output_bam + " -o " + output_sam_sort;
        run_command(sam_sort_cmd, true);

        std::string sam_depth_cmd = "samtools depth " + output_sam_sort + " > " +
                                    output_sam_depth;
        run_command(sam_depth_cmd, true);

        std::string filter_FP_log = fp_dir + "/FP.log";
        std::string filter_FP_output = fp_dir + "/FP_reconstructe_seq";

        write_filtered_tr(output_sam_depth, reconstructed_seq_path,
                          filter_FP_output, filter_FP_log);

        //std::string save_org_cmd("mv " + reconstructed_seq_path + " " +
        //                         fp_dir + "/reconstructed_org.fasta");
        //run_command(save_org_cmd, true);
        //std::string rename_new_seq_cmd("mv " + filter_FP_output + " " +
                                       reconstructed_seq_path);
        //run_command(rename_new_seq_cmd, true);

        exit(0);

    }
    else
    {
        std::cout << "input does not contain pair or filter_paired parameter "
                  << "in config file set to false. So skip filter false positive "
                  << "with Hisat" << std::endl;
    }
    */
}

void Multi_graph_handler::
write_filtered_tr(std::string depth_file, std::string in_tr_file,
                                    std::string out_tr_file, std::string log_file)
{
    typedef tsl::hopscotch_map<std::string, uint32_t> Tr_hits_map;
    typedef tsl::hopscotch_map<std::string, uint32_t>::iterator Tr_hits_map_iterator;
    Tr_hits_map tr_hits;

    std::ifstream depth_file_reader(depth_file);
    std::string field0_str, field1_str;
    while(  std::getline(depth_file_reader, field0_str, '\t') &&
            std::getline(depth_file_reader, field1_str))
    {
        Tr_hits_map_iterator it = tr_hits.find(field0_str);
        if(it != tr_hits.end())
        {
            it.value()++;
        }
        else
        {
            tr_hits.insert(std::make_pair(field0_str, 1));
        }
    }
    depth_file_reader.close();

    //typedef tsl::hopscotch_map<uint32_t, uint32_t> Tr_len_map;
    //typedef tsl::hopscotch_map<std::string, uint32_t>::iterator Tr_len_map_iterator;
    //Tr_len_map tr_len;


    std::ofstream file_writer(out_tr_file);
    std::string header, read_base;

    std::ifstream read_file_reader(in_tr_file);
    while(  std::getline(read_file_reader, header) &&
            std::getline(read_file_reader, read_base))
    {
        //tr_len.insert(std::make_pair(header, read_base.size()));
        std::string seq_name = header.substr(1);

        Tr_hits_map_iterator it = tr_hits.find(seq_name);
        uint32_t hit_num = 0;
        if(it != tr_hits.end())
            hit_num = it->second;
        else
            hit_num = 0;
        if(hit_num >= read_base.size() * pair_fp_thresh)
        {
            file_writer << header << std::endl;
            file_writer << read_base << std::endl;
        }
    }
}

void Multi_graph_handler::reverse_complement(std::string & seq)
{
    std::reverse(seq.begin(), seq.end());
    for(std::string::iterator it=seq.begin(); it!=seq.end(); it++)
    {
        if(*it=='A')
            *it = 'T';
        else if(*it=='T')
            *it = 'A';
        else if(*it=='C')
            *it = 'G';
        else if(*it=='G')
            *it = 'C';
        else
        {
            shc_log_error("find unknown char %c\n", *it);
            exit(1);
        }
    }
}
