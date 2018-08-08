#include "Contig_handler.h"

bool contig_handler_sort (contig_num_t i,contig_num_t j) { return (i<j); }

struct Block_timer ch_timer;

void Contig_handler::dump_all_contig(std::string & filename)
{
    shc_log_info(shc_logname, "Start Dumping all contigs into file %s\n",
                                           filename.c_str());
    int reminder_bases_num = 0;
    std::ofstream outfile (filename.c_str());
    char four_base[4];
    //first line is about the size of contig
    if(is_use_compress)
        outfile << "c\t";
    else
        outfile << "n\t";

    outfile << contig_list.size() << std::endl;
    int j = 0;

    {
        std::string message = "Dumping " + std::to_string(num_contig)
                               +  " contig into file\n";
        Progress_bar progress{std::cout, 70u, message};
        uint64_t progress_step = num_contig/100;

        uint64_t num_contig_processed = 0;

        for(int i=0; i<num_contig; i++)
        {
            outfile << ">Contig_"<< (j++) << "\t" << (int)mean_count[i] << std::endl;
            // -1 since exclusive at the end
            //shc_log_info(shc_logname, "delimitor.at(i) : %d\n",delimitor.at(i));
            //shc_log_info(shc_logname, "delimitor.at(i+1) : %d\n",delimitor.at(i+1));

            for(int j=delimitor.at(i); j<delimitor.at(i+1); j++)
            {
                if(is_use_compress)
                {
                    if(contig_len_list[i]%FOUR_BASE==0 || j<delimitor.at(i+1)-1)
                    {
                        decode(contig_list.at(j), FOUR_BASE, four_base);
                        for(int k=0; k<FOUR_BASE; k++)
                            outfile << four_base[k];
                    }
                    else
                    {
                        reminder_bases_num = contig_len_list[i]%FOUR_BASE;
                        decode(contig_list.at(j), reminder_bases_num, four_base);
                        for(int k=0; k<reminder_bases_num; k++)
                            outfile << four_base[k];
                    }
               }
               else
               {
                    outfile << contig_list[j];
               }
            }
            outfile << std::endl;

            num_contig_processed++;
            if(num_contig_processed%progress_step ==0)
            {
                double percentage = static_cast<double>(num_contig_processed) /
                                    static_cast<double>(num_contig);
                progress.write(percentage);
            }
        }
    }
    shc_log_info(shc_logname, "Finish Dumping all contigs into file\n");
    outfile.close();
}

void Contig_handler::load_contig_file(std::string & filename)
{
    shc_log_info(shc_logname, "Reading contigs from the file\n");
    std::ifstream fileReader (filename.c_str());
    std::string line, temp, mc_s, is_comp, list_size;
    int byte_list_size = 0;
    size_t contig_compress_len;
    kmer_count_t mc;

    //read byte size
    std::getline(fileReader, is_comp,'\t') && std::getline(fileReader, list_size);

    byte_list_size = std::stoi(list_size);
    if(byte_list_size !=0 )
    {
        contig_list.reserve(byte_list_size);
        contig_list.resize(byte_list_size);
    }
    else
    {
        shc_log_error("contig file contains 0 contig\n");
    }


    start_timer(&ch_timer);
    while (std::getline(fileReader, temp, '\t') &&
            std::getline(fileReader, mc_s))
    {
        mc = std::stoi(mc_s);

        std::getline(fileReader, line);
        //std::cout << " lc "<<line <<std::endl;
        //after this, line must contain the real seq
        if(line[0] == '>')
        {
            shc_log_error("contig file wrong format, seq starts with >\n");
            exit(1);
        }

        if(is_use_compress)
        {
            contig_compress_len = encode_base_string(line.c_str(),
                           &(contig_list.at(delimitor.back())), line.size());
            delimitor.push_back(delimitor.back()+contig_compress_len);
        }
        else
        {
            memcpy(&(contig_list.at(delimitor.back())),
                                                   line.c_str(), line.size());
            delimitor.push_back(delimitor.back()+ line.size());
        }
        num_contig++;
        mean_count.push_back(mc);
        contig_len_list.push_back(line.size());
    }
    //print_delimitor();
    //print_contig_length();
    //print_contig(0);

    //std::cout << "Finish contig loading, ";
    //stop_timer(&ch_timer);

    fileReader.close();
    shc_log_info(shc_logname, "Finish reading contig %u from the file\n",
                                                num_contig);
}
/**
 * Accepting this contig and update information
 * @param mean_c
 * @param contig_len
 */
void Contig_handler::declare_new_contig(kmer_count_t mean_c, size_type contig_len)
{
    size_t byte_lenght = 0;
    if(is_use_compress)
    {
        uint8_t * contig_start = &(contig_list.at(delimitor.back()));
        byte_lenght = encode_in_place((char*)contig_start, contig_start, contig_len);
        //note the last is empty, and is reserved for delimitor to start
        contig_list.resize(delimitor.back()+byte_lenght);
        num_contig++;
        delimitor.push_back(delimitor.back()+byte_lenght);
        mean_count.push_back(mean_c);
        contig_len_list.push_back(contig_len);
    }
    else
    {
        contig_list.resize(delimitor.back() + contig_len);
        if(MAX_CONTIG_NUM == num_contig)
        {
            shc_log_error("type currently used is 32bit, is unable to hold number of contigs"
                    "change the type contig_num_t and Macro MAX_CONTIG_NUM to uint64_t and UINT64_MAX "
                    "under src/shc_type.h\n");
            exit(0);
        }

        num_contig++;
        delimitor.push_back(delimitor.back()+ contig_len);
        mean_count.push_back(mean_c);
        contig_len_list.push_back(contig_len);
    }

#ifdef LOG_CONTIG
    shc_log_info(shc_logname, "contig number increase to %u\n", num_contig-1);
    shc_log_info(shc_logname, "delimitor adds an interval %u->%u \n",
                               delimitor[delimitor.size()-2], delimitor[delimitor.size()-1]-1);
    shc_log_info(shc_logname, "new contig has length %d\n", contig_len);
    shc_log_info(shc_logname, "new contig has byte length %d\n", byte_lenght);
    shc_log_info(shc_logname, "*****************************\n");
#endif
}

/**
 * In case, contig is pushed in the list, but is not needed, resize the list
 * back to previous length.
 */
void Contig_handler::reject_new_contig()
{
    contig_list.resize(delimitor.back());
    //shc_log_info(shc_logname, "reject contig, last delimitor remains %u\n", delimitor.back());
}

void Contig_handler::print_contig_length()
{
    for(size_type i=0; i< num_contig; i++)
        std::cout << contig_len_list.at(i) << " ";
    std::cout << std::endl;
}

void Contig_handler::print_contig(contig_num_t i)
{
    for(size_type j=delimitor.at(i); j<delimitor.at(i+1); j++)
        std::cout << (uint16_t)contig_list[j] << " ";
    std::cout << std::endl;
}

void Contig_handler::log_contig(contig_num_t i, size_t len)
{
    uint8_t * contig_start = &(contig_list.at(delimitor[i]));
    *(contig_start + len) = '\0';
    shc_log_info(shc_logname, "contig: %s \n", (char*) contig_start);
}

void Contig_handler::print_delimitor()
{
    for(std::vector<size_type>::iterator i = delimitor.begin(); i!=delimitor.end(); i++)
    {
        std::cout << *i << ' ';
    }
    std::cout << std::endl;
}


void Contig_handler::update_delimitor_and_mean_count(contig_num_t * remove_list, contig_num_t remove_num)
{
    //construct the list
    std::vector<size_type> length_list;
    std::vector<size_type> new_length_list;
    std::vector<size_type> new_delimitor;
    std::vector<kmer_count_t> new_mean_count;

    for(size_type i=0; i<num_contig; i++)
    {
        length_list.push_back(delimitor[i+1]-delimitor[i]);
    }

    // remove the contig in length array
    contig_num_t remove_counter = 0;
    contig_num_t coming_remove_contig = remove_list[remove_counter];
    for(contig_num_t i=0; i<num_contig; i++)
    {
        //std::cout << "coming_remove_contig is " << coming_remove_contig<< std::endl;
        if (i!=coming_remove_contig)
        {
            new_length_list.push_back(length_list[i]);
            new_mean_count.push_back(mean_count[i]);
        }
        else
        {
            ++remove_counter;
            if(remove_counter<remove_num)
                coming_remove_contig = remove_list[remove_counter];
            else
                coming_remove_contig = num_contig;
        }
    }

    new_delimitor.push_back(0);
    for(std::vector<size_type>::iterator i=new_length_list.begin();
            i<new_length_list.end(); i++)
    {
        new_delimitor.push_back(new_delimitor.back()+ *i );
    }
    delimitor.swap(new_delimitor);
    delimitor.resize(num_contig+1-remove_num);
    delimitor.shrink_to_fit();

    mean_count.swap(new_mean_count);
    mean_count.resize(num_contig-remove_num);
    mean_count.shrink_to_fit();

    //free memory
    std::vector<size_type>().swap(new_delimitor);
    std::vector<size_type>().swap(new_length_list);
    std::vector<size_type>().swap(length_list);
    std::vector<kmer_count_t>().swap(new_mean_count);
}
