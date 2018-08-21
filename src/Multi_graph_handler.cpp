#include "Multi_graph_handler.h"

int Multi_graph_handler::
count_num_component_sparse_flow()
{
    num_comp = 0;
    for(auto& entry : boost::make_iterator_range(
        boost::filesystem::directory_iterator(setting.local_files.output_seq_graph_path), {}))
    {
        int test_index = setting.local_files.output_seq_graph_path.size() + 7;
        std::string entry_str = entry.path().string();
        //std::cout << "entry_str " << entry_str << std::endl;
        if(entry_str[test_index] == 'p')
            num_comp++;
        //std::cout << entry << "\n";
    }
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
    std::string cmd("cat " + all_files_path + " > " + setting.local_files.reconstructed_single_path);
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
    std::string cmd("cat ");
    for(int i=0; i<num_parallel; i++)
    {
        std::string file_path =
                    setting.local_files.output_seq_graph_result_path +
                    "/process"+std::to_string(i);
        cmd += file_path + " ";
    }
    cmd += " > " + setting.local_files.reconstructed_sf_path;
    //std::cout << "cmd " << cmd << std::endl;
    run_command(cmd, true);

    //rm all temp process
    for(int i=0; i<num_parallel; i++)
    {
        std::string rm_cmd("rm ");
        std::string file_path =
                    setting.local_files.output_seq_graph_result_path +
                    "/process"+std::to_string(i);
        rm_cmd += file_path;
        run_command(rm_cmd, false);
    }
    std::cout << "remove all temp process sparse flow output file under comp_graph_output" << std::endl;
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
    for(int i=0; i<work_list.size(); ++i)
    {
        size_t file_size;
        if(get_read_file_size(i, file_size))
            work_size_list.emplace_back(work_list.at(i), file_size);
        else
        {
            shc_log_error("%d file not exists\n", i);
            exit(1);
        }
    }

    File_sorter file_sorter;
    std::sort(work_size_list.begin(), work_size_list.end(), file_sorter);
    for(int i=0; i<work_size_list.size(); i++)
    {
        work_list[i] = work_size_list[i].first;
    }
}

void Multi_graph_handler::
partition_work_to_process(int num_thread, std::deque<int> & work_list,
                                std::vector<std::deque<int>> & process_queue)
{
    std::vector<int> place_order;
    for(int i=0; i<num_thread; i++ )
        place_order.push_back(i);
    for(int i=num_thread-1; i>=0; i-- )
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

bool Multi_graph_handler::get_read_file_size(int i, size_t & size)
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
    std::deque<int> work_queue;
    std::deque<int> * work_queue_ptr = &(work_queue);
    std::vector<std::deque<int>> process_queue(num_parallel, std::deque<int>());

    int start_class = 0;
    int end_class = num_components;
    for(int i=start_class; i<end_class; i++)
    {
        work_queue.push_back(i);
    }
    sort_work_list_by_size(work_queue);
    partition_work_to_process(num_parallel, work_queue, process_queue);
    add_or_overwrite_directory(setting.local_files.output_seq_graph_path,
                               setting.local_files.output_path);
    if(!setting.local_files.single_node_dir.empty())
    {
        //std::cout <<"single_node_dir " << setting.local_files.single_node_dir << std::endl;
        add_or_overwrite_directory(setting.local_files.single_node_dir,
                                       setting.local_files.output_path);
    }
    for(int i=0; i<process_queue.size(); i++)
    {
        //std::cout << "process " << i << ": ";
        std::deque<int> & my_works = process_queue[i];
        std::cout << "process " << i << " has " << my_works.size()
                  << " components to process" << std::endl;;
        info_log_info(shc_logname, "process %d: ", i);
        for(std::deque<int>::iterator it=my_works.begin();
                            it!=my_works.end(); it++)
        {
        //    std::cout << (*it) << " ";
            info_log_info(shc_logname, "%d ", (*it));
        }
        //std::cout << std::endl;
        info_log_info(shc_logname, "\n");
    }

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
        {
            usleep(5000);
            continue;
        }


        std::cout << "Child " << getpid() << std::endl; // Child can keep going and fork once

        // children process starts here
        int process_queue_index = i;
        std::deque<int> my_works = process_queue[process_queue_index];
        auto engine = std::default_random_engine{};
        std::shuffle(std::begin(my_works), std::end(my_works), engine);
        pthread_mutex_t work_lock;
        if (pthread_mutex_init(&work_lock, NULL) != 0)
            { shc_log_error("unable to initialize mutex\n"); exit(1);}


        Seq_graph_works seq_graph_work;
        seq_graph_work.work_list = &my_works;
        seq_graph_work.setting = setting;
        seq_graph_work.work_lock_ptr = &work_lock;
        seq_graph_work.init_total_work = my_works.size();
        seq_graph_work.process_i = process_queue_index;

        int num_threads = 2;
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
           // else
           //{
           //   std::cout << "created thread " << i << std::endl;
           //}
        }

        for(int i=0; i<threads.size(); i++)
        {
           pthread_join(threads[i], NULL);
        }
        pthread_mutex_destroy(&work_lock);
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

        int num_threads = 2;
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
    Block_timer timer;
    start_timer(&timer);

    std::vector<std::string> headers;
    std::vector<std::string> seqs;

    std::string header, seq, line, rmer;
    rmer.resize(rmer_length);
    uint64_t contig_num = 0;
    uint64_t byte;
    std::ifstream file_reader(in_file);
    uint64_t num_seq;

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

    rmer_contig_map.reserve(num_seq*ave_read_length);

    while( std::getline(file_reader, line))
    {
        if( line[0] == '>')
        {
            header = line.substr(1);
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
            if(line.size() > setting.output_seq_min_len)
            {
                seqs.push_back(line);
                for(uint64_t i=0; i<line.size()-rmer_length+1; i++)
                {
                    encode_kmer(&line.at(i), &byte, rmer_length);
                    Seq_info seq_info(contig_num, i);
                    Rmer_contig_map_iterator it = rmer_contig_map.find(byte);
                    if(it != rmer_contig_map.end())
                    {
                        (it.value()).push_back(seq_info);
                    }
                    else
                    {
                        std::vector<Seq_info> vec_seq(1, seq_info);
                        vec_seq.reserve(VEC_INIT_SIZE);
                        rmer_contig_map.insert(std::make_pair(byte, vec_seq));
                    }
                }
                contig_num++;
            }
            else
            {
                headers.pop_back();
            }
        }
    }
    file_reader.close();

    contig_num = 0;
    std::ofstream file_writer(output_file);
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

    std::vector<Seq_info> & first_contig_vec = first_rc_it.value();
    std::vector<Seq_info> & last_contig_vec = last_rc_it.value();

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
            shc_log_error("find unknown char\n");
            exit(1);
        }
    }
}
