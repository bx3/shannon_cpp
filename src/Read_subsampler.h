#ifndef READ_SUBSAMPLER_H
#define	READ_SUBSAMPLER_H

#include <stdlib.h>
#include <time.h>
#include "shc_type.h"
#include "shc_google_sparsehash.h"

#define SIMPLE_SAMPLE 0
#define PROB_SAMPLE 1
#define SAMPLE_SPACE 1000000

struct Read_subsampler {
    Read_subsampler( bool is_store_all_reads_and_features_,
        int read_sampler_k_, unsigned int seed )
                    : read_sampler_k(read_sampler_k_),
                      k(0), is_always_pass(false), subsample_factor(0),
                      is_store_all_reads_and_features(is_store_all_reads_and_features_)
    {
        if(read_sampler_k_ == 0)
            is_always_pass = true;
        else if(read_sampler_k < 0)
        {
            subsample_factor = -read_sampler_k;
            is_use_prob_strategy = false;
        }
        else
        {
            k = read_sampler_k;
            is_use_prob_strategy = true;
        }
        if(seed == 0 )
        {
            srand(time(NULL));
        }
        else
        {
            srand(seed);
            std::cout << "read sampler set seed " << seed << std::endl;
        }
    }

    void setup(int num_comp, unsigned int seed)
    {
        num_read_seen.assign(num_comp, 0);
        num_read_accepted.assign(num_comp, 0);
    }

    inline bool simple_decider(const int & comp_i)
    {
        if(is_always_pass)
        {
            return true;
        }

        if(subsample_factor <= 0)
        {
            shc_log_error("simpler decider current subsample_factor %d, check for bug\n",
                                                        subsample_factor);
            exit(1);
        }

        if(num_read_seen[comp_i]++ % subsample_factor == 0)
        {
            num_read_accepted[comp_i]++;
            //shc_log_info(shc_logname, "seen %d, accepted %d, subsample_factor %d\n",
            //                    num_read_seen, num_read_accepted, subsample_factor);
            return true;
        }
        return false;
    }

    inline bool decide_to_keep_read(const int & comp_i,
                            kmer_count_t * kmer_counts, int num_kmer)
    {
        if(is_always_pass)
        {
            return true;
        }

        if(is_use_prob_strategy)
        {
            double avg_count = 0;
            for(int i=0; i<num_kmer; i++)
                avg_count += kmer_counts[i];
            avg_count /= num_kmer;

            if (avg_count == 0)
                return false;
            double prob_to_sample_read = k/avg_count;
            if(prob_to_sample_read >= 1.0)
                return true;
            else
                return (rand() % SAMPLE_SPACE < SAMPLE_SPACE*prob_to_sample_read);
        }
        else
        {
            return  simple_decider(comp_i);
        }
    }

    //return true if want to keep comp
    inline bool decide_to_keep_read(double & avg_count, const int & comp_i, double & prob_to_sample_read)
    {
        if(is_always_pass)
        {
            prob_to_sample_read = 1.0;
            return true;
        }

        if(is_use_prob_strategy)
        {
            if (avg_count == 0)
                return false;

            prob_to_sample_read = k/avg_count;
            if(prob_to_sample_read >= 1.0)
                return true;
            else
                return (rand() % SAMPLE_SPACE < SAMPLE_SPACE*prob_to_sample_read);
        }
        else
        {
            return  simple_decider(comp_i);
        }
    }

    double k;
    bool is_always_pass;
    bool is_use_prob_strategy;

    bool is_store_all_reads_and_features;

    int read_sampler_k;
    int subsample_factor;

    std::vector<uint64_t> num_read_seen;
    std::vector<uint64_t> num_read_accepted;
};

#endif
