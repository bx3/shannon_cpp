#include "json_parser.h"

std::string get_setting_string(Shannon_C_setting & setting)
{
    std::string setting_str;

    setting_str += "kmer_length is           : " + std::to_string(setting.kmer_length) + "\n";
    setting_str += "output_seq_min_len is    : " + std::to_string(setting.output_seq_min_len) + "\n";
    return setting_str;
}

int64_t format_memory_arg(std::string input)
{
    std::string units = "kMGT";
    char unit = input.back();
    std::size_t found = units.find(unit);
    if(found == std::string::npos)
    {
        for(int i=0; i<input.size(); i++)
            assert(input[i] >= '0' && input[i] <= '9');
        return boost::lexical_cast<int64_t>(input);
    }
    else
    {
        for(int i=0; i<input.size()-1; i++)
            assert(input[i] >= '0' && input[i] <= '9');
        std::string num_str = input.substr(0,input.size()-1);
        int64_t num = boost::lexical_cast<int64_t>(num_str);
        if(unit == 'k')
            return num * 1024;
        else if (unit == 'M')
            return num * 1024*1024;
        else if (unit == 'G')
            return num * 1024*1024*1024;
        else if (unit == 'T')
            return num * 1024*1024*1024*1024;
        else
        {
            std::cout << "unknown type " << input << std::endl;
            exit(1);
        }
    }
}

void parser_setting_file(std::string & file_path, Shannon_C_setting & setting)
{
    using boost::property_tree::ptree;

    std::ifstream json_file(file_path);
    ptree pt;
    boost::property_tree::json_parser::read_json(json_file, pt);
    // top level

    uint8_t kmer_length = pt.get<uint8_t>("kmer_length");
    bool is_double_stranded = pt.get<bool>("double_strand");
    bool has_pair = pt.get<bool>("has_pair");
    bool has_single = pt.get<bool>("has_single");
    bool is_compress = pt.get<bool>("is_compress");
    read_length_t single_read_length = pt.get<read_length_t>("single_read_length");
    read_length_t pair_1_read_length = pt.get<read_length_t>("pair_1_read_length");
    read_length_t pair_2_read_length = pt.get<read_length_t>("pair_2_read_length");
    //bool filter_single = pt.get<bool>("filter_single");
    //bool filter_paired = pt.get<bool>("filter_paired");
    int output_seq_min_len = pt.get<int>("output_seq_min_len");
    bool take_single_node_seq = pt.get<bool>("take_single_node_seq");
    bool take_contig_seq = pt.get<bool>("take_contig_seq");


    //int num_parallel = pt.get<int>("num_parallel");

    setting.kmer_length = kmer_length;
    setting.has_single = has_single;
    setting.has_pair = has_pair;
    setting.is_double_stranded = is_double_stranded;
    setting.is_compress = is_compress;
    setting.single_read_length = single_read_length;
    setting.pair_1_read_length = pair_1_read_length;
    setting.pair_2_read_length = pair_2_read_length;
    //setting.filter_single = filter_single;
    //setting.filter_paired = filter_paired;
    setting.output_seq_min_len = output_seq_min_len;
    setting.take_single_node_seq = take_single_node_seq;
    setting.take_contig_seq = take_contig_seq;

    ptree multi_graph_setup = pt.get_child("multi_graph_setup");
    setting.num_parallel = multi_graph_setup.get<int>("num_parallel");

    std::string avail_mem_str = multi_graph_setup.get<std::string>("avail_mem");
    int64_t avail_mem =  format_memory_arg(avail_mem_str);
    setting.avail_mem = avail_mem;

    ptree find_rep = pt.get_child("find_rep");
    setting.rmer_length = find_rep.get<int>("rmer_length");

    // file system
    ptree file_structure = pt.get_child("file_structure");

    std::string input_jf_path   = file_structure.get<std::string>("input_jf_path");
    std::string input_kmer_path = file_structure.get<std::string>("input_kmer_path");
    // OK if both empty or both filled
    if(input_kmer_path.empty() && !input_jf_path.empty())
    {
        shc_log_error("kmer dict and jf files must be provided together\n");
        exit(1);
    }
    else if(!input_kmer_path.empty() && input_jf_path.empty())
    {
        shc_log_error("kmer dict and jf files must be provided together\n");
        exit(1);
    }

    // read
    std::string input_read_path = file_structure.get<std::string>("input_read_path");
    std::string input_read_path_1 = file_structure.get<std::string>("input_read_path_1");
    std::string input_read_path_2 = file_structure.get<std::string>("input_read_path_2");
    // output dir
    std::string output_path = file_structure.get<std::string>("output_path");
    std::string reference_seq_path = file_structure.get<std::string>("reference_seq_path");

    boost::filesystem::path base_path_boost = boost::filesystem::current_path();
    std::string base_path(base_path_boost.string());

    Local_files lf(has_single, has_pair,
                   base_path, output_path,
                   input_kmer_path, input_read_path,
                   input_read_path_1, input_read_path_2, reference_seq_path,
                   input_jf_path);

    setting.local_files = lf;


    // dup_correction
    ptree dup_correction = pt.get_child("kmer_setup");
    uint8_t rmer_length = dup_correction.get<uint8_t>("rmer_length");
    double threshold = dup_correction.get<double>("threshold");
    kmer_count_t min_count = dup_correction.get<kmer_count_t>("min_count");
    size_t min_length = dup_correction.get<size_t>("min_length");
    bool is_use_set = dup_correction.get<bool>("is_use_set");
    double load_factor = dup_correction.get<double>("load_factor");
    int num_sort_thread = dup_correction.get<int>("num_sort_thread");
    setting.dup_setting.set_para(rmer_length, threshold, min_count, min_length,
        is_use_set, load_factor, num_sort_thread);

    // contig setup
    ptree contig_graph_setup = pt.get_child("contig_graph_setup");
    int num_test = contig_graph_setup.get<idx_t>("num_test");
    bool is_assign_best = contig_graph_setup.get<bool>("is_assign_best");
    int read_sampler_k = contig_graph_setup.get<int>("read_sampler_k");
    setting.contig_graph_setup.set_para(num_test, is_assign_best, read_sampler_k);

    //metis_setup
    ptree metis_setup = pt.get_child("metis_setup");
    bool use_multiple_partition = metis_setup.get<bool>
                                                    ("use_multiple_partition");
    idx_t partition_size = metis_setup.get<idx_t>("partition_size");
    idx_t non_partition_size = metis_setup.get<idx_t>("non_partition_size");
    idx_t penalty = metis_setup.get<idx_t>("penalty");
    real_t overload = metis_setup.get<real_t>("overload");
    setting.metis_setup.set_para(use_multiple_partition, partition_size,
                                 overload,  penalty, non_partition_size);
    // seq_graph_setup
    ptree seq_graph_setup = pt.get_child("seq_graph_setup");
    int hop_pair = 0;//seq_graph_setup.get<int>("max_hop_for_pair_search"); //not used anymore
    int hop_path = seq_graph_setup.get<int>("max_hop_for_known_path");
    int mate_pair_len = 300;//seq_graph_setup.get<int>("mate_pair_len");  //not used anymore
    setting.seq_graph_setup.set_para(hop_pair, hop_path, mate_pair_len);

    //sparse_flow setup
    ptree sparse_flow_setup = pt.get_child("sparse_flow_setup");
    int multiple_test = sparse_flow_setup.get<int>("multiple_test");
    setting.sparse_flow_setup.set_para(multiple_test);
    shc_log_info(shc_logname, "parsed json info\n");
}
