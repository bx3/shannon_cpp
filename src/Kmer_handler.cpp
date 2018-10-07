    #include "Kmer_handler.h"

size_t kmer_curr_proc = 0;
size_t kmer_prev_proc = 0;

Block_timer kh_timer;

/**
 * Constructor
 * @param length : kmer length
 */

void Kmer_handler::run_kmer_handler()
{
    start_timer(&kh_timer);
    std::cout << "start sorting Kmer using count" << std::endl;
    sort_kmer_descending_count();
    log_stop_timer(&kh_timer);
#ifdef SHOW_PROGRESS
    std::cout << "finish sorting Kmer using count, ";
    stop_timer(&kh_timer);
#endif

    start_timer(&kh_timer);

    if(build_dict_from_kmer_file() != 0)
    {
        shc_log_error("Error reading kmers\n");
    }
    log_stop_timer(&kh_timer);

    //std::string first_level_filter_kmers = setting.local_files.output_path + "/first_level_filter_kmers";
    //dump_loaded_kmers_to_file(first_level_filter_kmers);

    start_timer(&kh_timer);
    find_contig();
    log_stop_timer(&kh_timer);

}

void Kmer_handler::load_kmer_into_dict()
{
    if(lf->input_kmer_path.empty())
    {
        lf->input_kmer_path = lf->algo_input + "/kmer.dict";
        if(!exist_path(lf->input_kmer_path))
        {
            std::cerr << ("Need to add kmer path to the argument or json file, it is empty now\n");
            exit(1);
        }
    }
    shc_log_info(shc_logname, "Start loads unfiltered kmer into dict at path %s\n",
                                lf->input_kmer_path.c_str());
    std::cout << "Start loads unfiltered kmer into dict at path: "
              << lf->input_kmer_path << std::endl;
    std::ifstream fileReader(lf->input_kmer_path);
    std::string kmer_base, count_s;
    uint64_t count;
    kmer_count_t encoded_count;
    uint64_t byte;

    Kmer_info kmer_info(0,0,false, IMPOSSIBLE_CONTIG_NUM);
    // for showing progress
    uint64_t num_kmer_loaded = 0;
    Block_timer read_timer;
    start_timer(&read_timer);
    std::cout << "num_dict_kmer " << num_dict_kmer << std::endl;
    RESERVE_NUM_KMER(kmer_counter, num_dict_kmer*1.5)

    std::cout << "load_kmer_into_dict reserve " << (num_dict_kmer*1.5)
              << " kmers" << std::endl;

    {
        std::string message = "Loading " + std::to_string(num_dict_kmer)
                            + " kmers into dict\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t progress_step = num_dict_kmer/100;

        start_timer(&kh_timer);
        while(  std::getline(fileReader, kmer_base, '\t') &&
                std::getline(fileReader, count_s))
        {
            if(num_kmer_loaded%progress_step ==0)
            {
                double percentage = static_cast<double>(num_kmer_loaded) /
                                    static_cast<double>(num_dict_kmer);
                progress.write(percentage);
                shc_log_info(shc_logname, "loaded %u kmers out of %u\n",
                                        num_kmer_loaded, num_dict_kmer);
                log_stop_timer(&read_timer);
                start_timer(&read_timer);
            }

            //std::cout << "kmer_base " << kmer_base << std::endl;
            //std::cout << "count_s " << count_s << std::endl;
            encode_kmer(kmer_base.c_str(), &byte, kmer_length);
            count = boost::lexical_cast<uint64_t>(count_s);
            encoded_count = encode_count(count, compress_ratio);

            if(encoded_count > KMER_MAX_COUNT)
            {
                shc_log_error("kmer count is too large \n");
                exit(1);
            }
            num_kmer_loaded++;

            kmer_info.count = encoded_count;
            kmer_counter[byte] = kmer_info;

            //if(setting.is_double_stranded)
            //{
            //    num_kmer_loaded++;
            //    encode_reverse_kmer(kmer_base.c_str(), &byte, kmer_length);
            //    complement_num(&byte, kmer_length);
            //    kmer_counter[byte] = kmer_info;
            //}
            //std::cout << kmer_base<< " " << contig <<" "<<count <<" "<< info << std::endl;
            //std::cout << kmer_base<< " " << kmer_counter[byte].contig <<" "
            //          <<kmer_counter[byte].count <<" "<< kmer_counter[byte].info << std::endl;
        }
    }
    std::cout.flush();
    usleep(COUT_BUFFER_FLUSH_SLEEP);
    std::cout << std::endl;
    usleep(COUT_BUFFER_FLUSH_SLEEP);
    init_kmer_size = kmer_counter.size();
#ifdef SHOW_PROGRESS
    std::cout << "Finish kmer without info loading, ";
    stop_timer(&kh_timer);
#endif

    fileReader.close();
    shc_log_info(shc_logname, "kmer_length %u\n", kmer_length);
    shc_log_info(shc_logname, "finish load kmer\n");
}

/**
 * write all kmer and count into file, separated by space
 * @param filename
 */
void Kmer_handler::dump_kmers_with_info_to_file(std::string & filename)
{
    shc_log_info(shc_logname, "Start dumps kmer\n");
    start_timer(&kh_timer);
    std::ofstream outfile (filename.c_str());
    //std::ofstream unused_outfile (filename+ "_unused");

    Kmer_counter_map_iterator it;

    char kmer_base[KMER_HOLDER_LEN];
    outfile << kmer_counter.size()-num_skip_kmer << std::endl;

    Block_timer dump_timer;
    start_timer(&dump_timer);

    {
        std::string message = "Dumping " + std::to_string(kmer_counter.size())
                               +  " kmer with info into file\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t progress_step = kmer_counter.size()/100;

        uint64_t num_kmer_processed = 0;
        for (it = kmer_counter.begin(); it != kmer_counter.end(); it++)
        {
            if(it->second.used)
            {
                decode_kmer(kmer_base, &(it->first), kmer_length);
                kmer_base[kmer_length] = '\0';
                outfile << kmer_base << "\t" << it->second.contig
                                     << "\t" << get_count(it->second.count, compress_ratio)
                                     << "\t" << (uint16_t)(it->second.info)  << std::endl;
            }
            //else
            //{
            //    unused_outfile <<kmer_base << "\t" << it->second.contig
            //                         << "\t" << get_count(it->second.count, compress_ratio)
            //                         << "\t" << (uint16_t)(it->second.info)  << std::endl;
            //}

            if(num_kmer_processed%progress_step ==0)
            {
                double percentage = static_cast<double>(num_kmer_processed) /
                                    static_cast<double>(kmer_counter.size());
                progress.write(percentage);

                shc_log_info(shc_logname, "processed %u kmers out of %u\n",
                                        num_kmer_processed, kmer_counter.size());
                log_stop_timer(&dump_timer);
                start_timer(&dump_timer);
            }
            num_kmer_processed++;
        }
    }
    std::cout.flush();
    usleep(COUT_BUFFER_FLUSH_SLEEP);
    std::cout << std::endl;
    usleep(COUT_BUFFER_FLUSH_SLEEP);
    outfile.close();
    shc_log_info(shc_logname, "Finish\n");
}

void Kmer_handler::dump_loaded_kmers_to_file(std::string & filename)
{
    Block_timer local_timer;
    shc_log_info(shc_logname, "Start dumps kmer\n");
    start_timer(&local_timer);
    std::ofstream outfile (filename.c_str());

    Kmer_counter_map_iterator it;

    char kmer_base[KMER_HOLDER_LEN];

    for (it = kmer_counter.begin(); it != kmer_counter.end(); it++)
    {
        decode_kmer(kmer_base, &(it->first), kmer_length);
        kmer_base[kmer_length] = '\0';
        outfile << kmer_base << "\t" << it->second.count
                             << std::endl;
    }
    outfile.close();
    shc_log_info(shc_logname, "Finish\n");
    std::cout << "finish dump plain kmer with count" << std::endl;
    stop_timer(&local_timer);
}


void Kmer_handler::load_kmers_with_info_from_file(std::string & filename)
{
    shc_log_info(shc_logname, "Start loads kmer\n");
    std::ifstream fileReader(filename.c_str());
    std::string kmer_base, contig_s, count_s, info_s, kmer_counter_size_s;
    uint64_t count;
    kmer_count_t encoded_count;
    uint8_t info;
    contig_num_t contig;
    uint64_t byte;
    Kmer_info kmer_info(0,0,true, IMPOSSIBLE_CONTIG_NUM);
    // for showing progress
    uint64_t num_kmer_loaded = 0;
    Block_timer read_timer;
    start_timer(&read_timer);

    std::getline(fileReader, kmer_counter_size_s);
    uint64_t kmer_counter_size = boost::lexical_cast<uint64_t>(kmer_counter_size_s);

    RESERVE_NUM_KMER(kmer_counter, kmer_counter_size*1.5)


    std::cout << "reserve " << (kmer_counter_size*1.5) << " kmers" << std::endl;

    {
        std::string message = "Loading " + std::to_string(kmer_counter_size) +  " kmer with info into dict\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t progress_step = kmer_counter_size/100;

        start_timer(&kh_timer);
        while(  std::getline(fileReader, kmer_base, '\t') &&
                std::getline(fileReader, contig_s, '\t') &&
                std::getline(fileReader, count_s, '\t') &&
                std::getline(fileReader, info_s))
        {
            num_kmer_loaded++;
            if(num_kmer_loaded%progress_step ==0)
            {
                double percentage = static_cast<double>(num_kmer_loaded) /
                                    static_cast<double>(kmer_counter_size);
                progress.write(percentage);

                shc_log_info(shc_logname, "loaded %u kmers out of %u\n",
                                      num_kmer_loaded, kmer_counter_size);
                log_stop_timer(&read_timer);
                start_timer(&read_timer);
            }

            count = boost::lexical_cast<uint64_t>(count_s);
            encoded_count = encode_count(count, compress_ratio);
            if(encoded_count > KMER_MAX_COUNT)
            {
                shc_log_error("kmer count is too large \n");
                exit(1);
            }

            contig = boost::lexical_cast<contig_num_t>(contig_s);
            info = std::stoi(info_s);        // lexical_cast gives error

            encode_kmer(kmer_base.c_str(), &byte, kmer_length);
            kmer_info.contig = contig;
            kmer_info.info = info;
            kmer_info.count = encoded_count;
            kmer_counter[byte] = kmer_info;
            //std::cout << kmer_base<< " " << contig <<" "<<count <<" "<< info << std::endl;
            //std::cout << kmer_base<< " " << kmer_counter[byte].contig <<" "
            //          <<kmer_counter[byte].count <<" "<< kmer_counter[byte].info << std::endl;
        }
    }
    std::cout.flush();
    usleep(COUT_BUFFER_FLUSH_SLEEP);
    std::cout << std::endl;
    usleep(COUT_BUFFER_FLUSH_SLEEP);

//    std::cout << "Finish kmer with info loading, ";
//    stop_timer(&kh_timer);


    fileReader.close();
    shc_log_info(shc_logname, "kmer_length %u\n", kmer_length);
    shc_log_info(shc_logname, "finish load kmer\n");
}

void Kmer_handler::build_dict_from_kmer_file_helper(std::string file)
{
    std::string line_s, count_s;
    uint64_t count;
    uint64_t encoded_count;
    uint64_t kmer;
    Kmer_counter_map_iterator it;
    Kmer_info kmer_info;
    uint64_t num_kmer_loaded = 0;

    std::ifstream fileReader(file);
    //std::ofstream low_complex_writer(setting.local_files.output_path + "/low_complex_kmer");
    //std::cout << "Sorted kmer path " << lf->sorted_unfilter_file << std::endl;
    //shc_log_info(shc_logname, "load sorted kmer path %s\n", lf->sorted_unfilter_file.c_str());
    //std::cout << "sort file " << file << std::endl;

    Block_timer read_timer;
    start_timer(&read_timer);

    {
        std::string message = "Loading " + std::to_string(num_kmer) + " kmer into dict, discarding polyA kmer\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t progress_step = num_kmer/100;

        while (std::getline(fileReader, line_s, '\t') &&
           std::getline(fileReader, count_s) )
        {
            if(num_kmer_loaded%progress_step ==0)
            {
                double percentage = static_cast<double>(num_kmer_loaded) /
                                    static_cast<double>(num_kmer);
                progress.write(percentage);

                shc_log_info(shc_logname, "processed %u kmers out of %u\n",
                                        num_kmer_loaded, num_dict_kmer);
                log_stop_timer(&read_timer);
                start_timer(&read_timer);

            }

            if(is_polyA_del && is_low_complexity(line_s))
            {
                //std::cout << " low complex " << line_s << std::endl;
                //low_complex_writer << line_s << "\t"<< count_s << std::endl;
                continue;
            }

            count = boost::lexical_cast<uint64_t>(count_s);
            encoded_count = encode_count(count, compress_ratio);

            kmer_info.count = encoded_count;
            if (encoded_count > MAX_COUNT)
            {
                std::cout << "MAX_COUNT is  " << MAX_COUNT << std::endl;
                shc_log_error("kmer_count_t unable to hold count %d\n", encoded_count);
                exit(1);
            }
            encode_kmer(line_s.c_str(), &kmer, kmer_length);

            kmer_counter[kmer] = kmer_info;
            num_kmer_loaded++;
            if(setting.is_double_stranded)
            {
                num_kmer_loaded++;
                encode_reverse_kmer(line_s.c_str(), &kmer, kmer_length);
                complement_num(&kmer, kmer_length);
                kmer_counter[kmer] = kmer_info;
            }
        }
    }
    std::cout.flush();
    usleep(COUT_BUFFER_FLUSH_SLEEP);
    std::cout << std::endl;

    fileReader.close();
    //low_complex_writer.close();
#if defined(USE_DENSE_KMER) || defined(USE_SPARSE_KMER)
    kmer_counter.resize(kmer_counter.size());
#endif
}
/**
 * The function reads "kmer count" pair from file, each line denotes a
 * pair, and the delimitor is a tab
 * @param filename
 * @return
 */
int Kmer_handler::build_dict_from_kmer_file ()
{
    Block_timer local_timer;
    start_timer(&local_timer);

    shc_log_info(shc_logname, "Start load Kmer\n");
    build_dict_from_kmer_file_helper(lf->sorted_unfilter_file);
    shc_log_info(shc_logname, "finish kmer load %u\n", kmer_counter.size());

    init_kmer_size = kmer_counter.size();
    return 0;
}

// 1 is the first bit
bool Kmer_handler::is_info_ith_bit_set(uint8_t info, uint8_t i)
{
    if (i==1)
    {
        return (info & (uint8_t)(SHC_B1))>0;
    }
    else
    {
        return (info & (((uint8_t)(SHC_B1))<<(i-1)))>0;
    }
}

/**
 * This function checks the info in a kmer if it has prefix(1st) or suffix(2nd)
 * bit relating to another kmer.
 */
inline bool Kmer_handler::is_kmer_info_has_ps(uint8_t info)
{
    return info > 0;
}

/**
 *  This function checks the info in a kmer if it has prefix relating to another
 *  kmer.
 */
inline bool Kmer_handler::is_kmer_info_has_prefix(uint8_t info)
{
    return ((uint8_t)SHC_B1234 & info) > 0;
}

/**
 *  This function checks the info in a kmer if it has suffix relating to another
 *  kmer.
 */
inline bool Kmer_handler::is_kmer_info_has_suffix(uint8_t info)
{
    return ((uint8_t)SHC_B5678 & info) > 0;
}

void Kmer_handler::sort_kmer_descending_count_external()
{
    pid_t pid;

    if ((pid = fork()) < 0)
    {
            printf("fork error");
            _exit(1);
    } else if (pid == 0)
    {
            std::string parallel_arg("--parallel=");
            parallel_arg += std::to_string(dup_setting.num_sort_thread);

            if(dup_setting.sort_tmp_dir.empty())
            {
                std::string cmd(std::string("sort --parallel=") + std::to_string(dup_setting.num_sort_thread) +
                                        std::string(" -t '\t' -k2 -nr "));
                cmd += (lf->unfilter_file + " > " + lf->sorted_unfilter_file);
                print_yellow_cmd(cmd);
                if (execlp("sort", "sort", "-t","\t", "-k", "2", "-n", "-r",
                            parallel_arg.c_str(), lf->unfilter_file.c_str(),
                            "-o", lf->sorted_unfilter_file.c_str(),  (char *)0) < 0)
                {
                        printf("execlp error");
                        _exit(1);
                }
            }
            else
            {
                std::string cmd(std::string("sort --parallel=") + std::to_string(dup_setting.num_sort_thread) +
                                        std::string(" -t '\t' -k2 -nr ") + "-T " + dup_setting.sort_tmp_dir + " ");
                cmd += (lf->unfilter_file + " > " + lf->sorted_unfilter_file);
                print_yellow_cmd(cmd);
                if (execlp("sort", "sort", "-t","\t", "-k", "2", "-n", "-r",
                            parallel_arg.c_str(), "-T", dup_setting.sort_tmp_dir, lf->unfilter_file.c_str(),
                            "-o", lf->sorted_unfilter_file.c_str(),  (char *)0) < 0)
                {
                        printf("execlp error");
                        _exit(1);
                }
            }


            _exit(0);
    }

    pid_t profiler_pid;
    if( (profiler_pid = fork()) < 0)
    {
            printf("fork error");
            _exit(1);
    }
    else if (profiler_pid == 0)
    {
            printf("sort pid %d\n", pid);
            std::string pid_str = std::to_string(pid);

            if (execlp("./syrupy.py", "./syrupy.py", "-q", "-p", pid_str.c_str(), "-i", "10",
                        setting.local_files.mem_profiler.sort_log_path.c_str(),(char *)0) < 0)
            {
                    printf("execlp error");
                    _exit(1);
            }
            _exit(0);
    }

    if (waitpid(pid, NULL, 0) < 0)
    {
            printf("wait error");
            _exit(1);
    }
    if (waitpid(profiler_pid, NULL, 0) < 0)
    {
            printf("wait profiler error");
            _exit(1);
    }

}

/*
std::string num_thread_str= std::to_string(dup_setting.num_sort_thread);
std::string cmd(std::string("sort --parallel=") + num_thread_str + std::string(" -t '\t' -k2 -nr "));
cmd += (lf->unfilter_file + " > " + lf->sorted_unfilter_file);
run_command(cmd, true);
*/


/* sort by command
 sort -t$'\t' -k2 -nr small.dict > small_sort.dict
 */

void Kmer_handler::dump_and_sort_kmer_descending_count_external()
{
    shc_log_info(shc_logname,"start external sort\n");
    Block_timer local_timer;
    start_timer(&local_timer);


    shc_log_info(shc_logname,"after dump kmer\n");
    sort_kmer_descending_count_external();

    std::cout << "external sort finish " << std::endl;
    stop_timer(&local_timer);
}

/**
 * This function sorts kmer dictionary according to count information, and
 * store the descending list in a vector.
 */
void Kmer_handler::sort_kmer_descending_count()
{
    if (exist_path(lf->sorted_unfilter_file) )
    {
        if(get_filesize(lf->sorted_unfilter_file) > 0 )
        {
            std::string msg("Use existing sorted kmer file at path\n\t");
            msg += lf->sorted_unfilter_file + "\n";
            msg += "\tdelete this file if kmer.dict has changed.\n";
            print_important_notice(msg);
            return;
        }
    }

    lf->unfilter_file = lf->input_kmer_path;
    sort_kmer_descending_count_external();
}

/**
 * This function checks if the (K-1) suffix of a kmer has a common (k-1) prefix
 * with another kmer.
 * @param kmer_n      : the kmer of interest, if there is another kmer, it is
 *                      noted in the "info" field of the struct which can be
 *                      looked up by dictionary
 * @param new_kmer_n  : Similarly to above, that this another kmer has (k-1)
 *                      common prefix with another kmer, and is noted in "info".
 * @return            : 0 if this kmer has a suffix
 *                      -1 if this kmer does not have suffix
 */
bool Kmer_handler::find_suffix_kmer(uint64_t *kmer_n,
                                Kmer_counter_map_iterator & next_iter)
{
    uint64_t new_byte, best_byte;
    kmer_count_t count = 0;
    kmer_count_t max_count = 0;
    uint8_t best_base = DEFAULT_NUM;

    //write_all_suffix_info(kmer_n);
    Kmer_counter_map_iterator curr_kmer_iter = next_iter;
    write_all_prefix_info(kmer_n, curr_kmer_iter);
    Kmer_counter_map_iterator it_end = kmer_counter.end();

    for(uint8_t i=0; i<FOUR_BASE; i++)
    {
        new_byte = append_byte(kmer_n, i, kmer_length);

        Kmer_counter_map_iterator it = kmer_counter.find(new_byte);
        if (it != it_end)
        {
            write_kmer_info(i, false, curr_kmer_iter);   //record this kmer its possible prefix string
            if (!it->second.used)
            {
                count = it->second.count;
                /*
                if(count == max_count)
                {
                    std::cout << "find equal" << std::endl;
                    std::cout << "old char is " << num_to_char[best_base] << std::endl;
                    std::cout << "new char is " << num_to_char[i] << std::endl;
                }
                */
                if (count > max_count)
                {
                    best_base = i;
                    best_byte = new_byte;
                    max_count = count;
                    next_iter = it;
                    //shc_log_info(shc_logname, "max count %u with base %c\n",
                    //            max_count, num_to_char[best_base]);
                }
            }
        }
    }

    if (best_base == DEFAULT_NUM)
    {
        // so there is all kmer info entry is 0
        return false;
    }
    else
    {
        // if there are num_contig contigs, the last contig is numbered as num_contig-1
        *kmer_n = best_byte;
        ITER_SET_CONTIG(next_iter, ch->num_contig)

        return true;
    }
}

/**
 * This function checks if the (K-1) prefix of a kmer has a common (k-1) suffix
 * with another kmer.
 * @param kmer_n      : the kmer of interest, if there is another kmer, it is
 *                      noted in the "info" field of the struct which can be
 *                      looked up by dictionary
 * @param new_kmer_n  : Similarly to above, that if this another kmer has (k-1)
 *                      common suffix with the kmer of interest, then it is
 *                      noted in "info" field.
 * @return            : 0 if this kmer has a prefix
 *                      -1 if this kmer does not have prefix
 */
bool Kmer_handler::find_prefix_kmer(uint64_t *kmer_n,
                Kmer_counter_map_iterator & next_iter)
{
    char base_s[33];
    base_s[kmer_length] = '\0';
    decode_kmer(base_s, kmer_n, kmer_length);
    //shc_log_info(shc_logname, "input kmer %s\n", base_s);
    //char new_base_s[33];
    uint64_t new_byte, best_byte;
    kmer_count_t count = 0;
    kmer_count_t max_count = 0;
    uint8_t best_base = DEFAULT_NUM;

    //write_all_prefix_info(kmer_n);
    Kmer_counter_map_iterator curr_kmer_iter = next_iter;
    write_all_suffix_info(kmer_n, curr_kmer_iter);
    Kmer_counter_map_iterator it_end = kmer_counter.end();

    for(uint8_t i=0; i<FOUR_BASE; i++)
    {
        new_byte = prepend_byte(kmer_n, i, kmer_length);

        Kmer_counter_map_iterator it = kmer_counter.find(new_byte);

        if (it != it_end)
        {
            write_kmer_info(i, true, curr_kmer_iter);   //record this kmer its possible prefix string
            if(!it->second.used)
            {
                count = it->second.count;
                if (count > max_count)
                {
                    best_base = i;
                    best_byte = new_byte;
                    max_count = count;
                    next_iter = it;
                    //shc_log_info(shc_logname, "max count %u with base %c\n",
                    //            max_count, num_to_char[best_base]);
                }
            }
        }
    }

    if (best_base == DEFAULT_NUM)
    {
        // so there is all kmer info entry is 0
        //shc_log_info(shc_logname, "no output\n");
        return false;
    }
    else
    {
        *kmer_n = best_byte;
        ITER_SET_CONTIG(next_iter, ch->num_contig)
        return true;
    }
}

void Kmer_handler::write_all_suffix_info(const uint64_t *kmer_num,
                            Kmer_counter_map_iterator & curr_kmer_iter)
{
    uint64_t new_byte;
    Kmer_counter_map_iterator it;
    for(uint8_t i=0; i<FOUR_BASE; i++)
    {
        new_byte = append_byte(kmer_num, i, kmer_length);

        it = kmer_counter.find(new_byte);
        if(it != kmer_counter.end())
        {
            //record this kmer its possible prefix string
            write_kmer_info(i, false, curr_kmer_iter);
        }
    }
}

void Kmer_handler::write_all_prefix_info(const uint64_t *kmer_num,
                            Kmer_counter_map_iterator & curr_kmer_iter)
{
    uint64_t new_byte;
    Kmer_counter_map_iterator it;

    for(uint8_t i=0; i<FOUR_BASE; i++)
    {
        new_byte = prepend_byte(kmer_num, i, kmer_length);

        it = kmer_counter.find(new_byte);
        if(it != kmer_counter.end())
        {
            //record this kmer its possible prefix string
            write_kmer_info(i, true, curr_kmer_iter);
        }
    }
}

/**
 * This function is the only method which changes the info field
 * @param num           : one numeric number for one of "A T C G "
 * @param is_prefix     : is it a prefix, true for yes, false for no
 * @param kmer_num      : kmer in numeric form
 */
void Kmer_handler::write_kmer_info(uint8_t num, bool is_prefix,
                                Kmer_counter_map_iterator & kmer_it)
{
    //char base_string[KMER_HOLDER_LEN];
    //base_string[kmer_length] ='\0';
    //decode_kmer(base_string, kmer_num, kmer_length);
    uint8_t info = kmer_it->second.info;
    if (is_prefix)
    {
        info = info | (((uint8_t)SHC_B1)<<num);
    }
    else
    {
        info = info | ((uint8_t)SHC_B1<<(PREFIX_OFFSET+num));
    }

    ITER_SET_INFO(kmer_it, info)
}

/**
 * This function finds the contigs by going through all kmer in the dictionary.
 * For a specific kmer, it will find prefix first, then suffix. If it finds
 * a prefix, it will then start with the kmer which shares common (k-1)
 * character, and continue until there is no prefix kmer. Similarly for suffix.
 * At the end, it release all memory associated with sorting vector.
 * @return  number of contig
 */
contig_num_t Kmer_handler::find_contig()
{
    find_contig_external();
}

//return false if kmer is skipped
bool Kmer_handler::find_contig_helper(std::string & line_s, uint64_t count, bool is_rc)
{
    //shc_log_info(shc_logname, "find_contig_helper\n");
    char base_string[KMER_HOLDER_LEN];
    char root_string[KMER_HOLDER_LEN];
    char reverse_comp_root[KMER_HOLDER_LEN];
    base_string[kmer_length] = '\0';
    root_string[kmer_length] = '\0';
    reverse_comp_root[kmer_length] = '\0';
    uint64_t root_byte;// = boost::lexical_cast<uint64_t>(line_s);
    //char debug[1000000];
    if (is_rc)
    {
        encode_reverse_kmer(line_s.c_str(), &root_byte, kmer_length);
        complement_num(&root_byte, kmer_length);
        decode_kmer(reverse_comp_root, &root_byte, kmer_length);
    }
    else
    {
        encode_kmer(line_s.c_str(), &root_byte, kmer_length);
    }
    //std::cout << "line_s " << line_s << std::endl;
    //std::cout << "count_s " << count_s << std::endl;

    // start a new contig

    Kmer_counter_map_iterator root_iter = kmer_counter.find(root_byte);
    // needed since google dense hash still allow value access even key is deleted
    // updated needs command
    if(root_iter == kmer_counter.end())
        return true;

    if(count < dup_setting.min_count)
    {
        //shc_log_info(shc_logname, "preempt at iter %d, with count %u\n", i, GET_COUNT(root_iter->second.count, compress_ratio));
        //if(!root_iter->second.used)
        //    kmer_counter.erase(root_iter);
        num_skip_kmer++;
        return false;
    }


    Kmer_counter_map_iterator next_iter = root_iter;
    if (!root_iter->second.used)
    {
        ITER_SET_USED(root_iter, true)
        ITER_SET_CONTIG(root_iter, ch->num_contig)

        double contig_mean_count = 0.0;

        uint64_t total_count = count;//count;

        uint64_t total_kmer_in_contig = 1;
        uint64_t next_kmer = root_byte;
        //find prefix extension
        while (find_prefix_kmer(&next_kmer, next_iter))
        {

            ch->contig_list.push_back(decode_prefix_base(&next_kmer, kmer_length));
            total_count += get_count(next_iter->second.count, compress_ratio);

            total_kmer_in_contig++;

            ITER_SET_USED(next_iter, true)
        }
        // get everything in order
        std::reverse( ch->contig_list.begin() + ch->delimitor.back(),
                      ch->contig_list.end());

        // push the root string
        if(is_rc)
        {
            for(uint8_t t=0; t<kmer_length; t++)
                ch->contig_list.push_back(reverse_comp_root[t]);
        }
        else
        {
            for(uint8_t t=0; t<kmer_length; t++)
                ch->contig_list.push_back(line_s[t]);
        }
        //decode_kmer(root_string, &(root_byte), kmer_length);
        //shc_log_info(shc_logname, "Root kmer %s\n", line_s.c_str());
        //shc_log_info(shc_logname, "Root kmer %s\n", root_string);
        //shc_log_info(shc_logname, "kmer  : %s\n", line_s);
        next_kmer = root_byte;
        next_iter = root_iter;
        // find suffix extension
        //int suffix_len = 0;
        while (find_suffix_kmer(&next_kmer, next_iter) )
        {
            //place where info is recorded
            ch->contig_list.push_back(decode_suffix_base(&next_kmer));
            total_count += get_count(next_iter->second.count, compress_ratio);
            //shc_log_info(shc_logname, "kmer count %u\n",GET_COUNT(next_iter->second.count, compress_ratio));
            total_kmer_in_contig++;

            ITER_SET_USED(next_iter, true)
        }

        Contig_handler::size_type contig_len =
                            total_kmer_in_contig + kmer_length - 1;


        //shc_log_info(shc_logname, "total_count %d\n", total_count);
        contig_mean_count = static_cast<double>(total_count)/
                            static_cast<double>(total_kmer_in_contig);


        if(decide_contig_and_build_rmer(contig_mean_count, contig_len))
        {
            //uint8_t * contig_start = &(ch->contig_list.at(ch->delimitor[ch->num_contig]));
            //memcpy(debug, contig_start, contig_len);
            //debug[contig_len] = '\0';
            ch->declare_new_contig(std::ceil(contig_mean_count), contig_len);
            //std::cout << "accept, num contig " <<ch->num_contig << std::endl;
        }
        else
        {
            uint8_t * contig_start = &(ch->contig_list.at(ch->delimitor[ch->num_contig]));
            //restore_kmer_for_contig(contig_start, contig_len);
            delete_kmer_for_contig(contig_start, contig_len);
            ch->reject_new_contig();
            //std::cout << "reject, num contig " <<ch->num_contig << std::endl;
        }

        kmer_curr_proc += total_kmer_in_contig;
    }
    //else
    //{
    //    shc_log_info(shc_logname, "kmer used %u %s %d\n", ch->num_contig, line_s.c_str(), is_rc);
    //}
    return true;
}



contig_num_t Kmer_handler::find_contig_external()
{

    shc_log_info(shc_logname, "Start find contig\n");

    ch->contig_list.reserve(kmer_counter.size()/3);

    std::ifstream fileReader(lf->sorted_unfilter_file);

    std::string line_s, count_s;

    //int i = 0;
    start_timer(&prog_timer);

    {
        std::string message = "Connecting "+std::to_string(init_kmer_size) + " kmer into contig\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t progress_step = init_kmer_size/100;

        while (std::getline(fileReader, line_s, '\t') &&
           std::getline(fileReader, count_s) )
        {
            //std::cout << line_s << std::endl;
            //std::cout << count_s << std::endl;
            if(is_polyA_del && is_low_complexity(line_s))
            {
                continue;
            }

            uint64_t count = boost::lexical_cast<uint64_t>(count_s);
            if(!find_contig_helper(line_s, count, false))
                break;
            if(setting.is_double_stranded)
            {
                find_contig_helper(line_s, count, true);
            }

            if((kmer_curr_proc-kmer_prev_proc) > progress_step)
            {
                int percentage = (100 * kmer_curr_proc)/init_kmer_size;

                progress.write(static_cast<double>(kmer_curr_proc) /
                               static_cast<double>(init_kmer_size));

                log_stop_timer(&prog_timer);

                start_timer(&prog_timer);
                shc_log_info(shc_logname, "[ %u%] %u kmer out of %u is processed\n",
                        percentage, kmer_curr_proc, (init_kmer_size));
                kmer_prev_proc = kmer_curr_proc;
            }
        }
    }
    fileReader.close();
    std::cout.flush();
    usleep(COUT_BUFFER_FLUSH_SLEEP);
    std::cout << std::endl;
    usleep(COUT_BUFFER_FLUSH_SLEEP);

    if(dup_setting.is_use_set)
    {
        shc_log_info(shc_logname, "SET : rmer_set has %d entry, each with size %d\n",
                        rmer_set.size(), sizeof(uint64_t));
        std::cout << "SET : rmer_set size is " << rmer_set.size() << std::endl;
        //deallocate_rmer_set();
    }
    else
    {
        shc_log_info(shc_logname, "DICT: rmer_contig_map size is %d\n", rmer_contig_map.size());
        // for analysis
        uint64_t total_entry = 0;
        for(Rmer_contig_map_iterator co_it = rmer_contig_map.begin();
                                     co_it!= rmer_contig_map.end(); co_it++)
        {
            total_entry += co_it->second.size();
        }
        shc_log_info(shc_logname,"rmer_contig_map has %d entry, each with size %d\n",
                            total_entry, sizeof(contig_num_t));
        //deallocate_rmer_contig_map();
    }

    ch->contig_list.shrink_to_fit();
    //update the size, and release memory for kmer dict, see google sparsehash
#if defined(USE_DENSE_KMER) || defined(USE_SPARSE_KMER)
    kmer_counter.resize(0);
#endif
    shc_log_info(shc_logname, "kmer size is %llu\n", kmer_counter.size());
    shc_log_info(shc_logname, "Finish finding contig, %ld contigs, %d kmer_get deleted\n",
                                                   ch->num_contig, num_kmer_deleted);

    return ch->num_contig;
}

uint8_t Kmer_handler::num_bit_info(uint8_t info)
{
    uint8_t num = (uint8_t)(SHC_B1) & info;

    for(uint8_t i=1; i<8; i++)
    {
        num += (info >> i) & (uint8_t)(SHC_B1);
    }
    return num;
}

void Kmer_handler::traverse_kmer_count()
{
    char temp[33];
    temp[kmer_length] = '\0';
    shc_log_info(shc_logname, "Start traversing kmer counter dict\n");
    Kmer_counter_map_iterator it;
    uint8_t num_bit = 255;

    for (it = kmer_counter.begin(); it != kmer_counter.end(); it++)
    {
        num_bit = num_bit_info(it->second.info);
        //shc_log_info(shc_logname, "%u\n", num_bit);
        if(num_bit>2)
        {
            decode_kmer(temp, &(it->first), kmer_length);
            shc_log_info(shc_logname, "kmer %s in contig %u has more than two connectivity\n",
                                temp, it->second.contig);
        }
    }
    shc_log_info(shc_logname, "Finish traversing kmer counter dict\n");
}

/**
 * This function deallocates memory for sorting kmer vectors.
 */
void Kmer_handler::deallocate_kmer_descend_list()
{
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Start deallocate kmer sorting list\n");
#endif
    std::vector<Kmer_Occurence_Pair>().swap(kmer_descend_list);
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Finish deallocate kmer sorting list\n");
#endif
}

void Kmer_handler::deallocate_contig_count_map()
{
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Start deallocate contig count map for rmer list\n");
#endif
    ;
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Finish deallocate contig count map for rmer list\n");
#endif
}

void Kmer_handler::deallocate_rmer_contig_map()
{
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Start deallocate rmer contig map for rmer list\n");
#endif
    Rmer_contig_map().swap(rmer_contig_map);
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Finish deallocate rmer contig map for rmer list\n");
#endif
}

void Kmer_handler::deallocate_rmer_set()
{
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Start deallocate rmer set\n");
#endif
    Rmer_set().swap(rmer_set);
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Finish deallocate rmer set\n");
#endif
}

void Kmer_handler::clear_kmer_map()
{
    Kmer_counter_map().swap(kmer_counter);
}

void Kmer_handler::deallocate_kmer_map()
{
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Start deallocate kmer map\n");
#endif

    Kmer_counter_map().swap(kmer_counter);
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Finish deallocate kmer map\n");
#endif
}


void Kmer_handler::restore_kmer_for_contig(uint8_t * contig_start,
        Contig_handler::size_type contig_len)
{
    shc_log_info(shc_logname, "Start restore kmer\n");
    uint64_t byte = 0;
    char base_string[KMER_HOLDER_LEN];
    base_string[kmer_length] = '\0';
    Kmer_counter_map_iterator iter_end = kmer_counter.end();

    for(Contig_handler::size_type i=0; i<contig_len-kmer_length+1; i++)
    {
        encode_kmer((char*)(contig_start+i), &byte, kmer_length);
        Kmer_counter_map_iterator iter = kmer_counter.find(byte);
        if(iter == iter_end)
        {
            shc_log_error("see a unknown kmer\n");
            exit(1);
        }
        else
        {
            ITER_SET_USED(iter, false)
        }
    }
}

/**
 * The kmer are commented and is not considered used here for two reasons.
 * 1. Clearing the kmer map does not reduce the memory, and map does not
 * have a shrink_to_fit function.
 * 2. There are some problems to erase google sparsehash, it should have
 * a answer, but due to reason 1, it is not solved. When uncommenting, be aware.
 * @param i
 */
void Kmer_handler::delete_kmer_for_contig(uint8_t * contig_start,
        Contig_handler::size_type contig_len)
{
    //shc_log_info(shc_logname, "Start to delete kmer\n");
    uint64_t byte = 0;
    char base_string[KMER_HOLDER_LEN];
    base_string[kmer_length] = '\0';

    for(Contig_handler::size_type i=0; i<contig_len-kmer_length+1; i++)
    {
        encode_kmer((char*)(contig_start+i), &byte, kmer_length);
        kmer_counter.erase(byte);
    }
    //shc_log_info(shc_logname, "Finish delete kmer\n");
//#ifdef LOG_DELETED_KMER
    for(Contig_handler::size_type i=0; i<contig_len; i++)
    {
        info_log_info(lf->deleted_contig_path.c_str(), "%c", (char)contig_start[i]);
    }
    info_log_info(lf->deleted_contig_path.c_str(), "\n");
//#endif
    num_kmer_deleted += contig_len-kmer_length+1;

}

void Kmer_handler::get_contig_handler(Contig_handler * c_handler)
{
    ch = c_handler;
}

uint8_t Kmer_handler::get_kmer_length()
{
    return kmer_length;
}

/**
 * This function first checks if using a count or a list, then use contig length
 * and count to do a primary filter. If it passes, then check if it shares
 * common rmer with previous set. An exception is that if it is the first
 * rmer, we always accept the contig.
 * @param mean_count
 * @param contig len
 * @return if this contig is accepted, if not rmer is not built
 */
bool Kmer_handler::decide_contig_and_build_rmer(double mean_count, Contig_handler::size_type len)
{
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Start contig decision \t CANDIDATE CONTIG %d\n", ch->num_contig+1);
    shc_log_info(shc_logname, "Start length count filter\n");
#endif
    bool reject_flag =
         static_cast<double>(len) * std::pow(static_cast<double>(mean_count), 0.25) <
         2.0 * static_cast<double>(dup_setting.min_len) *
         std::pow(static_cast<double>(dup_setting.min_count), 0.25);
    if(len<dup_setting.min_len || reject_flag)
    {
#ifdef LOG_KMER
        shc_log_info(shc_logname, "len %d \n", len, mean_count);
        shc_log_info(shc_logname, "No pass, due to length %d count %d\n", len, mean_count);
#endif
        return false;
    }

    if(dup_setting.is_use_set)
    {
        return use_set_to_filter(mean_count, len);
    }
    else
    {
        return use_list_to_filter(mean_count, len);
    }
}

bool Kmer_handler::use_set_to_filter(kmer_count_t mean_count, Contig_handler::size_type len)
{
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Start set filter\n");
#endif
    uint8_t * contig_start = &(ch->contig_list.at(ch->delimitor[ch->num_contig]));
    char temp[32];
    temp[dup_setting.rmer_length] = '\0';
    uint64_t byte;

    uint64_t common_element = 0;
    Contig_handler::size_type total_rmer_num = len-dup_setting.rmer_length+1;

    //shc_log_info(shc_logname, "Calculating common element stat\n");
    std::vector<uint64_t> new_bytes;
    new_bytes.reserve(total_rmer_num);
    for(Contig_handler::size_type i=0; i<total_rmer_num; i++)
    {
        encode_kmer((char*)(contig_start+i), &byte, dup_setting.rmer_length);
        if(rmer_set.find(byte) != rmer_set.end())
            common_element++;
        else
            new_bytes.push_back(byte);
    }
#ifdef LOG_KMER
    shc_log_info(shc_logname, "\t\t\tSTAT summary\n");
    shc_log_info(shc_logname, "\t\t\tcontig length %d\n", len);
    shc_log_info(shc_logname, "\t\t\tCommon element %lld\n", common_element);
    shc_log_info(shc_logname, "\t\t\ttotal_rmer_num %lld\n", total_rmer_num);

    shc_log_info(shc_logname, "Start redundancy filter\n");
#endif
    if(static_cast<double>(common_element)/static_cast<double>(total_rmer_num)
                            < dup_setting.threshold)
    {
        for(std::vector<uint64_t>::iterator it= new_bytes.begin();
                                            it!=new_bytes.end(); it++)
        {
            rmer_set.insert(*it);
        }
#ifdef LOG_KMER
        shc_log_info(shc_logname, "Finish redundancy filter, Passing\n");
#endif
        return true;
    }
    else
    {
#ifdef LOG_KMER
        shc_log_info(shc_logname, "Finish redundancy filter, Fail to pass due to high redundance\n");
#endif
        return false;
    }
}

/**
 *
 * @param mean_count
 * @param len
 * @return  true when this contig is accepted, otherwise no
 */
bool Kmer_handler::use_list_to_filter(kmer_count_t mean_count, Contig_handler::size_type len)
{
#ifdef LOG_KMER
    shc_log_info(shc_logname, "Start list filter, with rmer %u\n", dup_setting.rmer_length);
#endif
    uint8_t * contig_start = &(ch->contig_list.at(ch->delimitor[ch->num_contig]));

    uint64_t byte;

    rmer_count_t max_till_now = 0;
    contig_num_t max_contig_index = IMPOSSIBLE_CONTIG_NUM;

    Contig_handler::size_type total_rmer_num = len-dup_setting.rmer_length+1;
    Contig_count_map contig_count_map;
    contig_count_map.set_empty_key(IMPOSSIBLE_CONTIG_NUM - 50);
    contig_count_map.set_deleted_key(IMPOSSIBLE_CONTIG_NUM - 2000);
    contig_count_map.resize(total_rmer_num);
    // finding the most similar previous contig
    for(Contig_handler::size_type i=0; i<total_rmer_num; i++)
    {
        //get rmer
        encode_kmer((char*)(contig_start+i), &byte, dup_setting.rmer_length);
        Rmer_contig_map_iterator rmer_contig_it = rmer_contig_map.find(byte);
        //if rmer corresponds to contig
        if(rmer_contig_it != rmer_contig_map.end())
        {
            std::vector<contig_num_t> & contig_vec = ITER_GET_VALUE(rmer_contig_it);
            //check the count and record in a dict
            for(std::vector<contig_num_t>::iterator it = contig_vec.begin();
                                               it != contig_vec.end(); it++)
            {
                contig_num_t contig_i = *it;
                Contig_count_map_iterator count_it = contig_count_map.find(contig_i);
                std::pair<Contig_count_map_iterator, bool> c_pair;
                rmer_count_t rmer_count = 1;
                if(count_it == contig_count_map.end())
                {
                    rmer_count = 1;
                    c_pair = contig_count_map.insert(std::make_pair(contig_i, rmer_count));
                }
                else
                {
                    count_it->second = count_it->second + 1;
                    rmer_count = count_it->second;
                    //rmer_count = ++(count_it->second);
                }
                /*
                if(rmer_count == max_till_now)
                {
                    shc_log_info(shc_logname, "equal count\n");
                }
                */
                if(rmer_count > max_till_now)
                {
                    max_till_now = rmer_count;
                    max_contig_index = contig_i; //which contig
                }
            }
        }
    }
#ifdef LOG_KMER
    shc_log_info(shc_logname, "rmer_contig_map size %d\n", rmer_contig_map.size());
    shc_log_info(shc_logname, "max contig is %d\n", max_contig_index);
#endif

    std::vector<uint8_t> indicator(len, 0);
    for(Contig_handler::size_type i=0; i<total_rmer_num; i++)
    {
        //get rmer
        encode_kmer((char*)(contig_start+i), &byte, dup_setting.rmer_length);
        Rmer_contig_map_iterator rmer_contig_it = rmer_contig_map.find(byte);
        //if rmer corresponds to contig
        if(rmer_contig_it != rmer_contig_map.end())
        {
            std::vector<contig_num_t> & contig_vec = ITER_GET_VALUE(rmer_contig_it);
            if(std::binary_search(contig_vec.begin(), contig_vec.end(), max_contig_index))
            {
                //shc_log_info(shc_logname, "binary search \n");
                //for(std::vector<contig_num_t>::iterator contig_vec_it=contig_vec.begin();
                //        contig_vec_it!=contig_vec.end(); ++contig_vec_it)
                //{
                //    shc_log_info(shc_logname, "%llu\n", *contig_vec_it);
                //}
                std::vector<uint8_t>::iterator rmer_start_it = indicator.begin()+i;
                std::fill(rmer_start_it, rmer_start_it+dup_setting.rmer_length, 0x01);
                //memset((void*)&(a.at(i)), 0x01, dup_setting.rmer_length);
            }
        }
    }

    int sum_of_elem = 0;
    for (uint8_t & elem : indicator)
        sum_of_elem += elem;

    if (static_cast<double>(sum_of_elem) > dup_setting.threshold * static_cast<double>(len))
    {
#ifdef LOG_KMER
        shc_log_info(shc_logname, "NO Pass, %d > %f\n", sum_of_elem,
                dup_setting.threshold * static_cast<double>(len));
#endif
        return false;
    }
    else
    {
#ifdef LOG_KMER
        shc_log_info(shc_logname, "Pass sum_of_elem: %d, len: %d\n", sum_of_elem, len);
#endif
        for(Contig_handler::size_type i=0; i<total_rmer_num; i++)
        {
            encode_kmer((char*)(contig_start+i), &byte, dup_setting.rmer_length);
            //shc_log_info(shc_logname, "dup_setting.rmer_length: %d, byte: %llx\n", dup_setting.rmer_length, byte);
            //since the current count has not been updated
            //also notice that if the accepted contig contains
            //the same kmer multiple times, the for the contig
            //list that corresponds to the kmer would contain
            //three element with the contig number
            //Change it if needed
            Rmer_contig_map_iterator rc_it = rmer_contig_map.find(byte);
            if(rc_it == rmer_contig_map.end())
            {
                rmer_contig_map.insert(std::make_pair(
                           byte, std::vector<contig_num_t>(1, ch->num_contig)));
            }
            else
            {
                if(rc_it->second.back()!= ch->num_contig)
                    ITER_GET_VALUE(rc_it).push_back(ch->num_contig);
            }
        }
        /*
        if(rmer_contig_map.empty())
        {
            shc_log_info(shc_logname, "contig_count_map empty\n");
        }
        else
        {
            for(Rmer_contig_map_iterator rc_it = rmer_contig_map.begin();
                                         rc_it != rmer_contig_map.end(); rc_it++ )
            {
                char base_r[32];
                base_r[15] = '\0';
                decode_kmer(base_r, &(rc_it->first), 15);
                shc_log_info(shc_logname, "base_r %s \n", base_r);
                for(std::vector<contig_num_t>::const_iterator v_it=rc_it->second.begin();
                        v_it != rc_it->second.end(); v_it++)
                {
                    shc_log_info(shc_logname, "contig %d \n", *v_it);
                }
            }
        }
        */
        //shc_log_info(shc_logname, "Pass\n");
        //shc_log_info(shc_logname, "Pass list filter, sum: %d, len: %d\n", sum_of_elem, len);
        return true;
    }
}

/**
* This function should only be used if kmer_length is set.
*/
uint8_t Kmer_handler::get_kmer_length_from_file(const std::string& filename)
{
    //std::assert(kmer_length != 0);
    std::string kmer_base, temp;
    std::ifstream file_reader(filename.c_str());
    std::getline(file_reader, temp); //the first line contain size info
    if(std::getline(file_reader, kmer_base, '\t') &&
       std::getline(file_reader, temp, '\t') &&
       std::getline(file_reader, temp, '\t') &&
       std::getline(file_reader, temp))
    {
        file_reader.close();
        return kmer_base.size();
    }
    else
    {
        shc_log_error("file is empty, or format wrong\n");
        exit(1);
    }
}


/**

    */
/*
 shc_log_info(shc_logname, "Start Filtering contigs\n");

    std::vector<contig_num_t> remove_list;
    for(contig_num_t i=0; i<num_contig; i++)
    {
        size_type len = delimitor[i+1] - delimitor[i];
        if(len<min_len || mean_count[i]<min_count)
        {
            shc_log_info(shc_logname, "Remove %u \n", i);
            remove_list.push_back(i);
        }
    }

    if(!remove_list.empty())
    {
        get_vector_string(remove_list, vec_str);
        shc_log_info(shc_logname, "Deleting vec %s\n", vec_str.str().c_str());
        delete_contig_list( &(remove_list.at(0)), remove_list.size());
    }
    else
    {
        shc_log_info(shc_logname, "length count checks and Nothing to delete\n");
    }

    remove_redundance(threshold);

    std::vector<contig_num_t>().swap(remove_list);
    shc_log_info(shc_logname, "Finish Filtering contigs\n");
 */

void Kmer_handler::load_sorted_kmer(std::string sorted_kmer_file)
{
    std::ifstream file_reader(sorted_kmer_file);

    std::string line_s, count_s;
    kmer_count_t count;
    uint64_t kmer;

    //int i = 0;

    while (std::getline(file_reader, line_s, '\t') &&
       std::getline(file_reader, count_s) )
    {
        count = boost::lexical_cast<kmer_count_t>(count_s);
        encode_kmer(line_s.c_str(), &kmer, kmer_length);
        kmer_descend_list.emplace_back(kmer, count);
        kmer_counter[kmer] = Kmer_info(count);
        //shc_log_info(shc_logname, "load sort kmer %s, %d\n", line_s.c_str(), count);
        //print_bit(kmer);

        //if(i++ > 10)
        //    exit(1);
    }
    file_reader.close();
    shc_log_info(shc_logname, "Finish load sorted kmer\n");

    /*
    for(Kmer_Occurence_Pair & item: kmer_descend_list )
    {
        Kmer_counter_map_iterator it = kmer_counter.find(item.first);
        if( it == kmer_counter.end())
        {
            shc_log_error("kmer counter does not have a kmer\n");
            exit(1);
        }
        else if(it->second.count != item.second)
        {
            shc_log_error("kmer counter and sorted list do not have same count\n");
            char base[33];
            base[kmer_length] = '\0';
            decode_kmer(base, &kmer, kmer_length);
            shc_log_error("%s \n", base);
            shc_log_error("kmer counter %d\n", it->second.count);
            shc_log_error("sorted list  %d\n", item.second);
            exit(1);
        }
    }
    */
}

bool Kmer_handler::is_low_complexity(std::string & base)
{
    int nA(0);
    int nC(0);
    int nG(0);
    int nT(0);
    for(int i=0; i<kmer_length; i++)
    {
        char letter = base.at(i);
        if (letter=='A')
        {
            nA++;
        }
        else if(letter=='C')
        {
            nC++;
        }
        else if(letter=='G')
        {
            nG++;
        }
        else if(letter == 'T')
        {
            nT++;
        }
    }
    int max_count = std::max(std::max(std::max(nA,nC),nG),nT);
    return max_count >= kmer_length - complex_thresh;
}

void Kmer_handler::log_kmer(uint64_t kmer)
{
    char base[33];
    base[kmer_length] = '\0';
    Kmer_info & kmer_info = kmer_counter[kmer];
    decode_kmer(base, &kmer, kmer_length);
    shc_log_info(shc_logname, "kmer %s\n", base);
    shc_log_info(shc_logname, "count %u, contig %u, info %u, used %s\n",
                            kmer_info.count, kmer_info.contig,
                            kmer_info.info, (kmer_info.used)?("Yes"):("No"));
    uint64_t next_byte;
    for(uint8_t i=0; i< BIT_PER_BYTE; i++)
    {
        if(IS_BIT_I_SET(kmer_info.info, i))
        {

            char next_base[33];
            next_base[kmer_length] = '\0';

            if(i<PREFIX_OFFSET)
            {
                next_byte = prepend_byte(&(kmer), i, kmer_length);
                decode_kmer(next_base, &next_byte, kmer_length);
                shc_log_info(shc_logname, "%u: %c prepended %s, %u, %s\n", i,
                                num_to_char[i],next_base, kmer_counter[next_byte].count,
                                (kmer_counter[next_byte].used)?("Yes"):("No"));
            }
            else
            {
                next_byte = append_byte(&(kmer), i-PREFIX_OFFSET, kmer_length);
                decode_kmer(next_base, &next_byte, kmer_length);
                shc_log_info(shc_logname, "%u: %c appended  %s, %u, %s\n", i,
                                num_to_char[i-PREFIX_OFFSET] ,next_base, kmer_counter[next_byte].count,
                                (kmer_counter[next_byte].used)?("Yes"):("No"));
            }
        }
    }
}

/*
 * int nA(0), nC(0), nG(0), nT(0);
    int n_remain(kmer_length);
    int thresh = kmer_length - complex_thresh;
    bool pA(true), pC(true), pG(true), pT(true);
    shc_log_info(shc_logname, "base %s\n", base.c_str()) ;
    for(int i=0; i<kmer_length; i++)
    {
        n_remain -= 1;
        char letter = base.at(i);
        if (letter=='A')
        {
            if(++nA + n_remain < thresh)
                pA = false;
        }
        else if(letter=='C')
        {
            if(++nC + n_remain < thresh)
                pC = false;
        }
        else if(letter=='G')
        {
            if(++nG + n_remain < thresh)
                pG = false;
        }
        else
        {
            if(++nT + n_remain < thresh)
                pT = false;
        }
        shc_log_info(shc_logname, "read %c\n", letter) ;
        shc_log_info(shc_logname, "pA %s\n", ((pA)?("YES"):("NO"))) ;
        shc_log_info(shc_logname, "pC %s\n", ((pC)?("YES"):("NO"))) ;
        shc_log_info(shc_logname, "pG %s\n", ((pG)?("YES"):("NO"))) ;
        shc_log_info(shc_logname, "pT %s\n", ((pT)?("YES"):("NO"))) ;
        if(!(pA || pC || pG || pT))
            return false;
    }
    return true;
 */
